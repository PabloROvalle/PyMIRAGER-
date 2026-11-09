"""Command-line interface and built-in presets."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from .io import read_case, write_outputs
from .forward import ForwardModel, parameters_for_case
from .inverse import RetrievalOptions, run_retrieval, corrlength


PRESETS = {'temperature': {'gaussian_width': 1.0, 'retrieval': {'filter_mode': 'none', 'c_temperature': 0.75, 'factor': 3.0, 'stop_mode': 'fortran'}}, 'abundance': {'gaussian_width': 0.925, 'retrieval': {'filter_mode': 'legacy_abundance', 'c_temperature': 3.0, 'c_molecule': 3.0, 'factor': 3.0, 'max_iterations': 20, 'rtol_chi2': 0.01, 'stop_mode': 'fortran'}}, 'abundance_split': {'gaussian_width': 0.925, 'retrieval': {'filter_mode': 'legacy_abundance', 'c_temperature': 0.75, 'c_molecule': 3.0, 'factor': 3.0, 'max_iterations': 20, 'rtol_chi2': 0.01, 'stop_mode': 'fortran'}}, 'aerosol': {'gaussian_width': 1.0, 'retrieval': {'filter_mode': 'legacy_aerosol', 'c_aerosol': 0.75, 'factor': 3.0}}, 'hcn': {'gaussian_width': 1.0, 'retrieval': {'filter_mode': 'legacy_hcn', 'c_special': 100.0, 'factor': 3.0}}, 'co2': {'gaussian_width': 1.2, 'retrieval': {'filter_mode': 'legacy_co2', 'c_special': 100.0, 'factor': 3.0}}, 'h2o': {'gaussian_width': 1.0, 'retrieval': {'filter_mode': 'legacy_h2o', 'c_special': 100.0, 'factor': 3.0}}, 'c6h6': {'gaussian_width': 0.925, 'retrieval': {'filter_mode': 'legacy_c6h6', 'c_special': 100.0, 'factor': 3.0}}, 'c6h6_fortran': {'gaussian_width': 0.925, 'retrieval': {'filter_mode': 'legacy_c6h6', 'c_temperature': 2.0, 'c_special': 2.0, 'factor': 3.0, 'rtol_chi2': 0.01, 'max_iterations': 20, 'stop_mode': 'fortran'}}, 'auto': {'gaussian_width': 0.925, 'retrieval': {'filter_mode': 'auto', 'factor': 3.0}}}

def options_for_preset(preset, overrides=None):
    config = dict(PRESETS[preset]['retrieval'])
    if overrides:
        config.update(overrides)
    return RetrievalOptions.from_dict(config)

def parser():
    p = argparse.ArgumentParser(description='Jupiter thermal IR radiative transfer and retrievals (5-20 um)')
    p.add_argument('--case', type=Path, required=True, help='Folder containing .fuel, .pta, .spe, resolution.txt and corrch4.txt')
    p.add_argument('--opacity-dir', type=Path, help='Folder containing jupiter_ch4_<gas>.opa')
    p.add_argument('--output', type=Path, default=Path('resultados_python'))
    p.add_argument('--preset', choices=PRESETS, default='auto', help='Retrieval preset (regularization, filters and instrumental response)')
    p.add_argument('--jacobian', choices=('original', 'finite_difference', 'analytical'), default='original', help='Molecular abundance Jacobian method (default: original)')
    p.add_argument('--options', type=Path, help='JSON file with optional retrieval settings')
    p.add_argument('--opacity-metadata', type=Path, help='JSON file defining alternate opacity grids')
    p.add_argument('--endian', choices=['little', 'big'], default='little')
    p.add_argument('--max-samples', type=int, default=160000)
    p.add_argument('--forward-only', action='store_true')
    p.add_argument('--check-inputs', action='store_true', help='Validate inputs without reading .opa files')
    p.add_argument('--without-ch4-legacy-scale', action='store_true', help='Disable the historical 7/6 scaling of the first gas')
    return p

def main(argv=None):
    a = parser().parse_args(argv)
    case = read_case(a.case, ch4_legacy_scale=not a.without_ch4_legacy_scale)
    print('Jupiter:', len(case.pressure), 'levels,', len(case.wave), 'spectral points')
    print('Gases:', ', '.join((f'{g.name}:{g.flag}' for g in case.gases)), '| cloudstate=', case.cloudstate)
    params = parameters_for_case(case)
    print('Retrieved parameters:', ', '.join((p.label for p in params)) or 'none')
    print('Molecular Jacobian:', a.jacobian)
    cfg = PRESETS[a.preset]
    overrides = json.loads(a.options.read_text()) if a.options else {}
    options = options_for_preset(a.preset, overrides)
    print(f'Preset {a.preset} | filters={options.filter_mode} | c(T)={options.c_temperature:g} | c(mol)={options.c_molecule:g} | factor={options.factor:g} | stopping={options.stop_mode}')
    if a.check_inputs:
        print('INPUT CHECK PASSED. No .opa files have been read or calculated.')
        return 0
    if not a.opacity_dir:
        a.opacity_dir = a.case
    model = ForwardModel(case, a.opacity_dir, case_dir=a.case, gaussian_width=cfg['gaussian_width'], max_samples=a.max_samples, opacity_metadata=a.opacity_metadata, endian='<' if a.endian == 'little' else '>', molecule_jacobian=a.jacobian)
    a.output.mkdir(parents=True, exist_ok=True)
    if a.forward_only:
        spec = model.forward()
        result = None
        chi2 = np.mean(((spec - case.observed) / case.error) ** 2)
        sigma = None
        ak = None
        dr = np.zeros(len(params))
        niter = 0
    else:
        def report_iteration(h):
            change = h['relative_chi2_change']
            change_text = 'N/A' if change is None else f'{change:+.6e}'
            print(f"iter {h['iteration']} mean chi2 {h['chi2_mean']:.6g} "
                  f"(chi2-old_chi2)/old_chi2 {change_text} "
                  f"max|delta| {h['max_abs_delta']:.5g}", flush=True)

        result = run_retrieval(model, options, progress=report_iteration)
        spec = result.synthetic
        chi2 = result.chi2
        sigma = result.sigma
        ak = result.averaging
        dr = result.dof
        niter = result.iterations
    write_outputs(case, params, spec, sigma, ak, dr, chi2, niter, a.output)
    if result is not None:
        (a.output / 'diagnostics.json').write_text(json.dumps({'chi2_mean': float(result.chi2), 'iterations': result.iterations, 'parameters': [p.label for p in params], 'correlation_length': [float(corrlength(p, options)) for p in params], 'filters': options.filter_mode, 'molecular_jacobian': a.jacobian, 'history': result.history}, indent=2, ensure_ascii=False))
    print(f'Results written to {a.output}; mean chi2 = {chi2:.6g}')
    return 0
