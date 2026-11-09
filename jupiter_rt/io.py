"""Input and output for Jupiter retrieval cases."""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from pathlib import Path
import re


@dataclass
class Gas:
    name: str
    scale: float
    flag: int

@dataclass
class Case:
    latitude: float
    mu_in: float
    mu_out: float
    cloudstate: int
    gases: list[Gas]
    pressure: np.ndarray
    temperature: np.ndarray
    profiles: np.ndarray
    taucloud: np.ndarray
    cloudmask: np.ndarray
    wave: np.ndarray
    observed: np.ndarray
    error: np.ndarray
    fwhm: np.ndarray

def read_fuel(file):
    raw = []
    for s in Path(file).read_text(errors='replace').splitlines():
        s = s.split('!')[0].split('#')[0].strip()
        if s:
            raw.append(s)
    if len(raw) < 5:
        raise ValueError('Incomplete .fuel file')
    lat, mu_in, mu_out = (float(raw[i].replace('D', 'E')) for i in range(3))
    cloudrow = raw[3].split()
    if cloudrow[0].lower() != 'cloud':
        raise ValueError('Expected cloud 0/1/2 in .fuel')
    cloud = int(cloudrow[1])
    gases = []
    for line in raw[4:]:
        fields = line.split()
        if len(fields) < 3:
            raise ValueError(f'Incomplete gas entry: {line}')
        gas = Gas(fields[0].lower(), float(fields[1].replace('D', 'E')), int(fields[2]))
        if gas.flag not in (0, 1, 2, 3):
            raise ValueError('Allowed gas flag: 0,1,2,3')
        gases.append(gas)
    if not gases:
        raise ValueError('No gas entries were provided')
    return (lat, mu_in, mu_out, cloud, gases)

def read_spe(file):
    with open(file, encoding='utf-8', errors='replace') as f:
        n = int(float(f.readline().split()[0]))
        dat = np.loadtxt(f)
    dat = np.atleast_2d(dat)
    if dat.shape != (n, 3):
        raise ValueError(f'.spe: expected {n} rows with 3 columns, got {dat.shape}')
    if np.any(dat[:, 2] <= 0):
        raise ValueError('Spectral uncertainties must be positive')
    return (dat[:, 0], dat[:, 1], dat[:, 2])

def read_resolution(file, nu):
    vals = np.loadtxt(file, skiprows=1)
    if vals.ndim != 2 or vals.shape[1] != 2:
        raise ValueError('resolution.txt: expected two columns')
    n = float(Path(file).read_text().splitlines()[0].strip())
    if int(round(n)) != len(vals):
        raise ValueError('resolution.txt row count does not match its header')
    x = vals[:, 0]
    y = vals[:, 1]
    o = np.argsort(x)
    if len(np.unique(x)) != len(x):
        raise ValueError('Repeated resolution wavenumbers')
    if min(nu) < min(x) or max(nu) > max(x):
        raise ValueError('Spectrum outside resolution.txt coverage; refusing to extrapolate')
    out = np.interp(nu, x[o], y[o])
    if np.any(out <= 0):
        raise ValueError('Instrumental resolution must be positive')
    return out

def read_case(path, *, ch4_legacy_scale=True):
    path = Path(path)
    lat, mui, mue, cloud, gases = read_fuel(path / 'jupiter_ch4.fuel')
    dat = np.loadtxt(path / 'jupiter_ch4.pta')
    # .pta layout: pressure, temperature, then one VMR column per .fuel gas.
    # For a special molecule, its VMR is the final column (no aerosol column).
    # Aerosol optical depth always starts from zero, independently of cloudstate.
    if dat.ndim != 2 or dat.shape[0] != 361:
        raise ValueError(f'.pta: expected 361 levels; got {dat.shape}')
    special_names = {'hcn', 'co2', 'h2o', 'c6h6'}
    if any(g.name in special_names for g in gases[:-1]):
        raise ValueError('.fuel: special molecule must be the last gas to match the final .pta VMR column')
    expected_columns = 2 + len(gases)
    has_special = gases[-1].name in special_names
    if dat.shape[1] == expected_columns:
        pass
    elif (not has_special and dat.shape[1] == expected_columns + 1
          and np.all(dat[:, -1] == 0)):
        # Accept a legacy unused trailing zero column, but never read it as aerosol.
        pass
    else:
        raise ValueError(
            f'.pta: expected {expected_columns} columns '
            f'(pressure, temperature and {len(gases)} gas VMRs); got {dat.shape[1]}. '
            'Aerosols are initialized to zero and do not have a .pta column.'
        )
    p = dat[:, 0].copy()
    temp = dat[:, 1].copy()
    profiles = dat[:, 2:2 + len(gases)].copy()
    if ch4_legacy_scale:
        profiles[:, 0] *= 7 / 6
    profiles *= np.array([g.scale for g in gases])[None, :]
    # Historical C6H6 abundance scaling from jupiter_inv_3D_c6h6/read_file.f90.
    # Keep the original input unchanged and apply the scale exactly once.
    for j, gas in enumerate(gases):
        if gas.name == 'c6h6':
            profiles[:, j] *= 50.0
    taucloud = np.zeros(len(p), dtype=float)
    mask = np.zeros(361)
    if cloud:
        fname = {1: 'nh3_clouds.txt', 2: 'ch4_hazes.txt'}.get(cloud)
        if fname is None:
            raise ValueError(f'Unsupported cloudstate: {cloud}')
        a = np.loadtxt(path / fname)
        if a.shape != (361, 2):
            raise ValueError('Cloud mask must contain 361 rows')
        if not np.allclose(p, a[:, 0], rtol=0.01):
            raise ValueError('Cloud-mask pressure grid does not match')
        mask = a[:, 1].copy()
    w, obs, err = read_spe(path / 'jupiter_ch4.spe')
    fw = read_resolution(path / 'resolution.txt', w)
    return Case(lat, mui, mue, cloud, gases, p, temp, profiles, taucloud, mask, w, obs, err, fw)

def write_outputs(case, params, synthetic, sigma, averaging, dr, chi2, niter, output_dir):
    from pathlib import Path
    path = Path(output_dir)
    path.mkdir(parents=True, exist_ok=True)
    np.savetxt(path / 'jupiter_ch4.res', np.column_stack((case.wave, case.observed, synthetic, case.error)), fmt='%.9e', header='nu_cm^-1 observed synthetic error')
    cols = [case.pressure, case.temperature]
    names = ['p_bar', 'T_K']
    for j, gas in enumerate(case.gases):
        if gas.flag in (2, 3):
            cols.append(case.profiles[:, j] / gas.scale)
            names.append(gas.name + '_profile_over_scale')
    if case.cloudstate:
        cols.append(case.taucloud)
        names.append('tau_cloud')
    if sigma is not None:
        for i, param in enumerate(params):
            cols.append(sigma[:, i])
            names.append('sigma_' + param.label)
    np.savetxt(path / 'jupiter_ch4.inv', np.column_stack(cols), fmt='%.9e', header=' '.join(names))
    if params and sigma is not None:
        np.savetxt(path / 'jupiter_ch4.sum', np.array([float(niter), float(chi2), *dr])[None, :], fmt='%.9e', header='iterations chi2_mean ' + ' '.join(('dof_' + x.label for x in params)))
        if averaging is not None:
            kernel = np.argmax(averaging, axis=0) + 1
            np.savetxt(path / 'jupiter_ch4.avg', np.column_stack([np.arange(1, 362), kernel]), fmt='%d', header='level ' + ' '.join((x.label for x in params)))
