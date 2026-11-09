"""Regularized inversion and retrieval iteration."""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass, field
from typing import Any
from scipy.linalg import cho_factor, cho_solve
from .forward import parameters_for_case


SPECIAL_MOLECULES = {'hcn', 'co2', 'h2o', 'c6h6'}

@dataclass
class RetrievalOptions:
    max_iterations: int = 15
    rtol_chi2: float = 0.01
    tol_step: float = 0.3
    stop_mode: str = 'combined'
    factor: float = 3.0
    factor_temperature: float | None = None
    factor_molecule: float | None = None
    factor_special: float | None = None
    factor_aerosol: float | None = None
    c_temperature: float = 0.75
    c_molecule: float = 3.0
    c_special: float = 100.0
    c_aerosol: float = 0.75
    aerosol_zero_from_index: int | None = None
    aerosol_zero_below_index: int | None = None
    filter_mode: str = 'auto'
    alpha_floor: float = 1e-26
    scale_hcn_delta: float = 0.1

    @classmethod
    def from_dict(cls, v):
        known = cls.__dataclass_fields__
        bad = set(v) - set(known)
        if bad:
            raise ValueError('Unknown options: ' + str(sorted(bad)))
        return cls(**v)

def filter_vectors(nlevel, width=15):
    idx = np.arange(nlevel, dtype=float)
    tropo = 1 - (np.tanh((idx - 90) / width) + 1) / 2
    strato = 1 - tropo
    return (tropo, strato)

def choose_filter_mode(params, mode):
    if mode != 'auto':
        return mode
    names = {p.label.lower() for p in params}
    if any((n in SPECIAL_MOLECULES for n in names)):
        if 'h2o' in names:
            return 'legacy_h2o'
        if 'hcn' in names:
            return 'legacy_hcn'
        if 'co2' in names:
            return 'legacy_co2'
        if 'c6h6' in names:
            return 'legacy_c6h6'
    if any((p.kind == 'aerosol' for p in params)):
        return 'legacy_aerosol'
    if any((p.kind == 'molecule' for p in params)):
        return 'legacy_abundance'
    return 'none'

def apply_filters(delta, params, options):
    delta = delta.copy()
    f, f2 = filter_vectors(delta.shape[0])
    mode = choose_filter_mode(params, options.filter_mode)
    valid = {'none', 'legacy_abundance', 'legacy_aerosol', 'legacy_hcn', 'legacy_co2', 'legacy_c6h6', 'legacy_h2o'}
    if mode not in valid:
        raise ValueError(f'Unknown filter: {mode}')
    if mode in ('legacy_abundance', 'legacy_co2'):
        if delta.shape[1] > 0:
            delta[:, 0] *= f
        if delta.shape[1] > 1:
            delta[:, 1] *= f2
    if mode in ('legacy_aerosol', 'legacy_c6h6'):
        if delta.shape[1] > 0:
            delta[:, 0] *= f
    if mode == 'legacy_hcn':
        if delta.shape[1] > 0:
            delta[:, 0] *= options.scale_hcn_delta
    return delta

def correlation(pressure, c):
    d = np.log(np.asarray(pressure))
    return np.exp(-0.5 * ((d[:, None] - d[None, :]) / c) ** 2)

def block_factor(param, options):
    if param.kind == 'temperature':
        return options.factor_temperature if options.factor_temperature is not None else options.factor
    if param.kind == 'aerosol':
        return options.factor_aerosol if options.factor_aerosol is not None else options.factor
    if param.label.lower() in SPECIAL_MOLECULES:
        return options.factor_special if options.factor_special is not None else options.factor
    return options.factor_molecule if options.factor_molecule is not None else options.factor

def corrlength(param, options):
    if param.kind == 'temperature':
        return options.c_temperature
    if param.kind == 'aerosol':
        return options.c_aerosol
    if param.label.lower() in SPECIAL_MOLECULES:
        return options.c_special
    return options.c_molecule

@dataclass
class StepResult:
    delta: np.ndarray
    sigma: np.ndarray
    averaging: np.ndarray
    dof: np.ndarray
    alpha: np.ndarray
    filters: str

def retrieval_step(pressure, kk, observed, synthetic, errors, params, options):
    nflux, nlevel, nparams = kk.shape
    if nparams != len(params):
        raise ValueError('Inconsistent number of Jacobian blocks')
    if nparams == 0:
        raise ValueError('No retrieval parameters selected')
    scovs = []
    alphas = []
    cov = np.diag(errors ** 2)
    trace_error = np.sum(errors ** 2) / np.sqrt(nflux)
    for j, param in enumerate(params):
        s = correlation(pressure, corrlength(param, options))
        k = kk[:, :, j]
        kt = k @ s
        trace_temp = np.sum(kt * k)
        if not np.isfinite(trace_temp) or trace_temp <= options.alpha_floor:
            raise ValueError(f'Jacobian {param.label}: has zero or degenerate sensitivity; check opacities and configuration')
        alpha = block_factor(param, options) * trace_error / trace_temp
        cov += alpha * kt @ k.T
        scovs.append(s)
        alphas.append(alpha)
    try:
        cf = cho_factor(cov, lower=True, check_finite=True)
        cy = cho_solve(cf, observed - synthetic)
        solve = lambda X: cho_solve(cf, X)
    except np.linalg.LinAlgError:
        cy = np.linalg.solve(cov, observed - synthetic)
        solve = lambda X: np.linalg.solve(cov, X)
    delta = np.zeros((nlevel, nparams))
    sigma = np.zeros_like(delta)
    A = np.zeros((nlevel, nlevel, nparams))
    dof = np.zeros(nparams)
    for j, param in enumerate(params):
        w = scovs[j] @ kk[:, :, j].T
        wi = solve(w.T).T
        delta[:, j] = alphas[j] * wi @ (observed - synthetic)
        A[:, :, j] = alphas[j] * wi @ kk[:, :, j]
        dof[j] = np.trace(A[:, :, j])
        sigma[:, j] = alphas[j] * np.sqrt(np.sum(wi ** 2 * errors[None, :] ** 2, axis=1))
    mode = choose_filter_mode(params, options.filter_mode)
    return StepResult(apply_filters(delta, params, options), sigma, A, dof, np.asarray(alphas), mode)

@dataclass
class RetrievalResult:
    synthetic: np.ndarray
    chi2: float
    iterations: int
    sigma: np.ndarray
    averaging: np.ndarray | None
    dof: np.ndarray
    history: list[dict[str, Any]] = field(default_factory=list)

def run_retrieval(model, options=None, *, progress=None):
    options = options or RetrievalOptions()
    case = model.case
    params = model.params
    for param in params:
        if param.kind == 'molecule' and np.any(case.profiles[:, param.gas_index] <= 0):
            raise ValueError(f'Profile {param.label} contains nonpositive abundances; cannot update ln(q).')
    if not params:
        spec = model.forward(with_jacobians=False)
        return RetrievalResult(spec, float(np.mean(((spec - case.observed) / case.error) ** 2)), 0, np.zeros((len(case.pressure), 0)), None, np.zeros(0))
    # Optional fixed-zero altitude boundaries for aerosol retrievals. These are
    # experimental constraints, not a physical truncation in the original Fortran.
    aerosol_lower = options.aerosol_zero_below_index
    aerosol_upper = options.aerosol_zero_from_index
    aerosol_columns = [j for j, param in enumerate(params) if param.kind == 'aerosol']
    aerosol_allowed = None
    if aerosol_lower is not None or aerosol_upper is not None:
        if not aerosol_columns:
            raise ValueError('Aerosol cutoffs require an aerosol retrieval')
        nlevel = len(case.pressure)
        if aerosol_lower is not None:
            if (isinstance(aerosol_lower, bool) or not isinstance(aerosol_lower, int)
                    or not 0 <= aerosol_lower < nlevel):
                raise ValueError('aerosol_zero_below_index must be an integer in [0, nlevel - 1]')
        if aerosol_upper is not None:
            if (isinstance(aerosol_upper, bool) or not isinstance(aerosol_upper, int)
                    or not 1 <= aerosol_upper < nlevel):
                raise ValueError('aerosol_zero_from_index must be an integer in [1, nlevel - 1]')
        if aerosol_lower is not None and aerosol_upper is not None and aerosol_lower >= aerosol_upper:
            raise ValueError('aerosol_zero_below_index must be smaller than aerosol_zero_from_index')
        aerosol_allowed = np.ones(nlevel, dtype=bool)
        if aerosol_lower is not None:
            aerosol_allowed[:aerosol_lower] = False
        if aerosol_upper is not None:
            aerosol_allowed[aerosol_upper:] = False
        case.taucloud[~aerosol_allowed] = 0.0

    chi_old = np.inf
    prev_step = np.inf
    hist = []
    last = None
    for iteration in range(options.max_iterations + 1):
        syn, kk = model.forward(with_jacobians=True)
        if aerosol_allowed is not None:
            # Fixed-zero aerosol levels have no sensitivity in the Jacobian.
            for j in aerosol_columns:
                kk[:, ~aerosol_allowed, j] = 0.0
        chi = float(np.mean(((syn - case.observed) / case.error) ** 2))
        if not np.isfinite(chi):
            raise FloatingPointError('Non-finite chi2')
        if iteration == options.max_iterations:
            hist.append({'iteration': iteration, 'chi2_mean': chi, 'stopped': 'max_iterations'})
            break
        step = retrieval_step(case.pressure, kk, case.observed, syn, case.error, params, options)
        if aerosol_allowed is not None:
            for j in aerosol_columns:
                step.delta[~aerosol_allowed, j] = 0.0
                step.sigma[~aerosol_allowed, j] = 0.0
                step.averaging[~aerosol_allowed, :, j] = 0.0
        abs_delta = float(np.max(np.abs(step.delta)))
        # Fractional chi-squared improvement relative to the previous iteration.
        # Positive values indicate a decrease in chi-squared; undefined at iter 0.
        chi2_rel_change = ((chi_old - chi) / chi_old
                           if np.isfinite(chi_old) and chi_old != 0 else None)
        hist.append({'iteration': iteration, 'chi2_mean': chi,
                     'relative_chi2_change': chi2_rel_change,
                     'max_abs_delta': abs_delta, 'dof': step.dof.tolist(),
                     'filter': step.filters})
        if progress:
            progress(hist[-1])
        last = step
        rel = abs(chi - chi_old) / max(chi, 1e-300)
        if options.stop_mode == 'fortran':
            is_converged = rel < options.rtol_chi2
        elif options.stop_mode == 'combined':
            is_converged = rel < options.rtol_chi2 and abs_delta < prev_step and (abs_delta < options.tol_step)
        else:
            raise ValueError(f'Unknown stopping mode: {options.stop_mode}')
        if is_converged:
            hist[-1]['stopped'] = 'converged'
            break
        chi_old = chi
        prev_step = abs_delta
        if abs_delta > 10:
            step.delta *= 0.5
        for j, param in enumerate(params):
            incr = step.delta[:, j]
            if param.kind == 'temperature':
                case.temperature[:] += incr
            elif param.kind == 'molecule':
                case.profiles[:, param.gas_index] *= np.exp(np.clip(incr, -80, 80))
            elif param.kind == 'aerosol':
                case.taucloud[:] = np.maximum(0, case.taucloud + incr)
                if aerosol_allowed is not None:
                    case.taucloud[~aerosol_allowed] = 0.0
        if np.any(case.temperature <= 0):
            raise FloatingPointError('Retrieval generated T<=0 K; check step size and regularization')
    syn, k_final = model.forward(with_jacobians=True)
    if aerosol_allowed is not None:
        for j in aerosol_columns:
            k_final[:, ~aerosol_allowed, j] = 0.0
    chi = float(np.mean(((syn - case.observed) / case.error) ** 2))
    final_diagnostics = retrieval_step(case.pressure, k_final, case.observed, syn, case.error, params, options)
    if aerosol_allowed is not None:
        for j in aerosol_columns:
            final_diagnostics.sigma[~aerosol_allowed, j] = 0.0
            final_diagnostics.averaging[~aerosol_allowed, :, j] = 0.0
    sigma = final_diagnostics.sigma
    A = final_diagnostics.averaging
    dof = final_diagnostics.dof
    return RetrievalResult(syn, chi, len(hist) - 1, sigma, A, dof, hist)
