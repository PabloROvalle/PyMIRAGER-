"""Thermal emission forward model."""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from scipy.sparse import csr_matrix
from .constants import C
from .geometry import atmosphere_geometry
from .opacity import OpacitySet, load_temperature_table
from .continuum import calculate_cia
from .jacobians import _molecular_jacobian_exact


@dataclass(frozen=True)
class Parameter:
    kind: str
    label: str
    gas_index: int = -1

def parameters_for_case(case):
    result = []
    if any((g.flag in (1, 3) for g in case.gases)):
        result.append(Parameter('temperature', 'T'))
    result.extend((Parameter('molecule', g.name, j) for j, g in enumerate(case.gases) if g.flag in (2, 3)))
    if case.cloudstate:
        result.append(Parameter('aerosol', 'aerosol'))
    return result

def planck(nu, T):
    nu = np.asarray(nu, float)
    hck = C.planck * C.cvel_cm_s / C.boltzmann
    with np.errstate(over='ignore', under='ignore', invalid='raise'):
        return 2 * C.planck * C.cvel_cm_s ** 2 * nu ** 3 / np.expm1(nu * hck / T)

def gaussian_operator(nu, wave, fwhm, gaussian_width=0.925):
    nu = np.asarray(nu, float)
    rows = []
    cols = []
    vals = []
    for i, (v, fw) in enumerate(zip(wave, fwhm)):
        a = int(np.searchsorted(nu, v - 5 * fw, side='right'))
        b = int(np.searchsorted(nu, v + 5 * fw, side='left'))
        if a >= b:
            raise ValueError(f'No convolution samples near {v:g} cm-1')
        indexes = np.arange(a, b)
        w = np.exp2(-((nu[indexes] - v) / (fw / (2 * gaussian_width))) ** 2)
        s = w.sum()
        if s <= 0:
            raise ValueError(f'Convolution has zero weight near {v} cm-1')
        rows.extend([i] * len(indexes))
        cols.extend(indexes.tolist())
        vals.extend((w / s).tolist())
    return csr_matrix((vals, (rows, cols)), shape=(len(wave), len(nu)))

class ForwardModel:

    def __init__(self, case, opacity_dir=None, *, endian='<', opacity_metadata=None, gaussian_width=0.925, max_samples=160000, case_dir=None, opacity_set=None, cia_provider=None, cia_library=None, molecule_jacobian='original', fd_step=0.001):
        if molecule_jacobian not in ('original', 'finite_difference', 'analytical'):
            raise ValueError('Molecular Jacobian must be original, finite_difference or analytical')
        if fd_step <= 0 or not np.isfinite(fd_step):
            raise ValueError('fd_step must be finite and positive')
        self.molecule_jacobian = molecule_jacobian
        self.fd_step = float(fd_step)
        self.case = case
        self.case_dir = case_dir
        self.params = parameters_for_case(case)
        self.max_samples = max_samples
        self.gaussian_width = gaussian_width
        self.opacity = opacity_set if opacity_set is not None else OpacitySet([g.name for g in case.gases], opacity_dir, metadata=opacity_metadata, endian=endian)
        self.cia_provider = cia_provider
        self.cia_library = cia_library
        self.table = load_temperature_table()
        self.nu = self._internal_grid()
        self.convolution = gaussian_operator(self.nu, case.wave, case.fwhm, gaussian_width)
        self.nu_jac = self.nu[::5]
        self.convolution_jac = gaussian_operator(self.nu_jac, case.wave, case.fwhm, gaussian_width)
        if hasattr(self.opacity, 'validate_active_coverage'):
            self.opacity.validate_active_coverage([g.name for g in case.gases if g.flag in (2, 3)], case.wave)

    def _internal_grid(self):
        c = self.case
        lo = min(c.wave - 5 * c.fwhm)
        hi = max(c.wave + 5 * c.fwhm)
        if lo < 500 or hi > 2000:
            raise ValueError(f'Spectral window [{lo:.3f},{hi:.3f}] exceeds the shared grid 500-2000 cm-1')
        bw = 10.0
        step = 0.001
        first = int(np.floor((lo - 500) / bw))
        last = int(np.ceil((hi - 500) / bw))
        nblocks = last - first
        if nblocks < 1:
            nblocks = 1
        num = nblocks * 10000
        if num > self.max_samples:
            raise MemoryError(f'{num} samples: safety limit {self.max_samples}. Split the observed spectrum.')
        return 500 + (first * 10000 + np.arange(num)) * step

    def forward(self, with_jacobians=False):
        c = self.case
        nlevel = len(c.pressure)
        nu = self.nu
        z, g, mu = atmosphere_geometry(c.pressure, c.temperature, c.latitude, c.mu_out)
        cscale = C.avogadro * 100000.0 / (2.6868e+19 * 10000.0 * C.mean_molar_mass_kg_mol)
        delta_p = np.r_[c.pressure[:-1] - c.pressure[1:], c.pressure[-1]]
        columns = c.profiles * (delta_p * cscale / g)[:, None]
        if self.cia_provider is not None:
            tau = np.asarray(self.cia_provider(c.pressure, c.temperature, g, nu), dtype=float)
        else:
            import pathlib
            if self.case_dir is None:
                raise ValueError('case_dir is required to locate corrch4.txt')
            tau = calculate_cia(pathlib.Path(self.case_dir), c.pressure, c.temperature, g, nu, library=self.cia_library)
        if tau.shape != (nu.size, nlevel):
            raise ValueError('CIA must return an array of shape (nfreq,nlevel)')
        molpars = [p for p in self.params if p.kind == 'molecule']
        kappas = {p.label: np.zeros_like(tau) for p in molpars} if with_jacobians else {}
        for k in range(nlevel):
            for j, gas in enumerate(c.gases):
                op = self.opacity.lines(gas.name, k, nu, self.table[k], c.temperature[k])
                effect = op * columns[k, j]
                tau[:, k] += effect
                if with_jacobians and gas.flag in (2, 3):
                    kappas[gas.name][:, k] = effect
        tau[:, :-1] += c.taucloud[None, :-1]
        transmission = np.exp(-np.cumsum((tau / mu[None, :])[:, ::-1], axis=1)[:, ::-1])
        t_layer = (c.temperature[:-1] + c.temperature[1:]) / 2
        intensity = np.zeros(nu.size)
        for k in range(nlevel - 1):
            intensity += planck(nu, t_layer[k]) * (transmission[:, k + 1] - transmission[:, k])
        spec = np.asarray(self.convolution @ intensity)
        if not with_jacobians:
            return spec
        jac = np.zeros((len(c.wave), nlevel, len(self.params)))
        trans5 = transmission[::5, :]
        p5 = np.empty((nu[::5].size, nlevel - 1))
        for k in range(nlevel - 1):
            p5[:, k] = planck(self.nu_jac, t_layer[k])
        for ip, param in enumerate(self.params):
            if param.kind == 'temperature':
                gk = np.empty((len(c.wave), nlevel - 1))
                for k in range(nlevel - 1):
                    radiance = (planck(self.nu_jac, t_layer[k] + 1) - p5[:, k]) * (trans5[:, k + 1] - trans5[:, k])
                    gk[:, k] = self.convolution_jac @ radiance
                jac[:, 0, ip] = gk[:, 0] * 0.5
                jac[:, 1:-1, ip] = (gk[:, :-1] + gk[:, 1:]) * 0.5
                jac[:, -1, ip] = gk[:, -1] * 0.5
            elif param.kind == 'molecule' and self.molecule_jacobian != 'original':
                jac[:, :, ip] = _molecular_jacobian_exact(transmission, nu, t_layer, kappas[param.label], mu, self.convolution, method=self.molecule_jacobian, fd_step=self.fd_step)
            elif param.kind in ('molecule', 'aerosol'):
                accum = np.zeros(len(self.nu_jac))
                for k in range(nlevel - 2):
                    accum -= (p5[:, k] - p5[:, k + 1]) * trans5[:, k]
                    if param.kind == 'molecule':
                        multiplier = kappas[param.label][::5, k] / mu[k]
                    else:
                        multiplier = c.cloudmask[k] / mu[k]
                    jac[:, k, ip] = self.convolution_jac @ (accum * multiplier)
        jac[np.logical_or(jac > 1, jac < -1)] = 0
        return (spec, jac)
