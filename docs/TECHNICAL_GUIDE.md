# Jupiter RT Python — Technical Guide

**Documented version:** 0.4.4  
**Project status:** Operational modular Python implementation, with ongoing validation against the original Fortran programs.  
**Documentation cutoff:** November 9, 2026.

> This guide complements the repository's introductory `README.md`; it does not replace it. It describes features implemented in v0.4.4, rather than treating planned features as already available.

## Contents

1. [Purpose and scope](#1-purpose-and-scope)
2. [Code structure](#2-code-structure)
3. [Requirements and installation](#3-requirements-and-installation)
4. [Preparing a case](#4-preparing-a-case)
5. [Running the program](#5-running-the-program)
6. [Physical and numerical model](#6-physical-and-numerical-model)
7. [Retrieval and regularization](#7-retrieval-and-regularization)
8. [Available presets](#8-available-presets)
9. [Jacobians](#9-jacobians)
10. [Reading and interpreting results](#10-reading-and-interpreting-results)
11. [Troubleshooting](#11-troubleshooting)
12. [Validation status and remaining work](#12-validation-status-and-remaining-work)

See also the [input/output file format reference](FILE_FORMATS.md) and the [scientific validation report](VALIDATION.md).

## 1. Purpose and scope

`jupiter_rt` is a **pure-Python** implementation of a thermal-infrared emission model and atmospheric-profile retrieval system for Jupiter, derived from several original Fortran programs. It supports retrievals of temperature, molecular abundances, aerosols, and selected trace species, while preserving legacy configuration options and allowing alternative molecular Jacobian methods to be compared.

The general working range is **500–2000 cm⁻¹** (approximately 5–20 µm). The current implementation assumes **361 pressure levels** and has not been generalized to other vertical dimensions. Its internal spectral grid has a spacing of **0.001 cm⁻¹**, arranged in 10 cm⁻¹ blocks. Instrumental response is applied through Gaussian convolution using `resolution.txt`.

The repository does **not** contain binary `.opa` opacity tables. These must be supplied separately at runtime. No Fortran compilation is required: the collision-induced absorption (CIA) routines have also been translated into Python.

## 2. Code structure

```text
Jupiter_RT_Python/
├── jupiter_rt/
│   ├── __init__.py
│   ├── __main__.py      # python -m jupiter_rt entry point
│   ├── cli.py           # command-line arguments and presets
│   ├── constants.py     # physical and planetary constants
│   ├── geometry.py      # atmospheric geometry, altitude, gravity
│   ├── io.py            # reading input cases and writing outputs
│   ├── opacity.py       # binary opacity tables and spectral grids
│   ├── cia_python.py    # H2-H2, H2-He, and H2-CH4 CIA components
│   ├── continuum.py     # continuum optical depth
│   ├── forward.py       # thermal emission, convolution, Jacobians
│   ├── jacobians.py     # alternative molecular Jacobians
│   ├── inverse.py       # retrieval iterations and regularization
│   └── data/
│       └── temperature_table.npy  # required bundled lookup table
├── docs/
│   ├── TECHNICAL_GUIDE.md
│   ├── FILE_FORMATS.md
│   └── VALIDATION.md
└── pyproject.toml
```

**Changing the target molecule does not require editing Python source**, provided that `.fuel` selects the correct species and that the corresponding `.opa` table is available. A preset specifies regularization, vertical filtering, and instrumental convolution; **it does not change the species or retrieval flags in `.fuel`**.

## 3. Requirements and installation

- Linux (the platform used for validation) with **Python 3.10 or later**.
- Declared dependencies: `numpy>=1.24`, `scipy>=1.10`.
- The relevant case files and `.opa` tables.
- Sufficient memory for the spectral grid and Jacobian arrays.

Install the package in the same Conda environment used to run it:

```bash
conda activate jwst_neptune
cd ~/MIRAGER/Jupiter_RT_Python
python3 -m pip install -e .
python3 -c "import jupiter_rt; print(jupiter_rt.__file__)"
python3 -m jupiter_rt --help
```

`~/MIRAGER/Jupiter_RT_Python` is the location used on the validation server, **not a required installation path**. After editable installation (`pip install -e .`), the module may be invoked from the case directory. The package also provides the `jupiter-rt` console command, but the examples use `python3 -m jupiter_rt` to make the Python interpreter explicit.

## 4. Preparing a case

A case consists of these files:

```text
case_directory/
├── jupiter_ch4.fuel
├── jupiter_ch4.pta
├── jupiter_ch4.spe
├── resolution.txt
└── corrch4.txt
```

For aerosol modes `cloud 1` or `cloud 2`, the case additionally requires `nh3_clouds.txt` or `ch4_hazes.txt`, respectively. Opacity tables may live in a different directory, specified using `--opacity-dir`.

The input reader expects these exact filenames. Column definitions, flags, and scaling conventions are described in [FILE_FORMATS.md](FILE_FORMATS.md).

**Fundamental `.pta` convention:** column 1 is pressure (bar), column 2 is temperature (K), and there is **one VMR column per gas, in the same order as `.fuel`**. There is no aerosol column: `taucloud` starts as an array of zeros. For HCN, CO₂, H₂O, and C₆H₆, the special molecule must be the **last gas** in both `.fuel` and `.pta`.

For backward compatibility, v0.4.4 also accepts one extra *all-zero* trailing column **only when no special molecule is present**; that column is ignored.

### Opacity tables

`--opacity-dir` must contain one `jupiter_ch4_<gas>.opa` file for every gas listed in `.fuel`, including gases that are not being retrieved. Names are matched literally: **`c2h3` remains the opacity identifier even when the species is scientifically referred to as C₂H₂**. It is not renamed automatically.

Validation-server example:

```text
/obs/tfouchet/data/opacite/jupiter/jwst
```

Do not use obsolete `/poubelle/...` paths that may appear in old Fortran logs as an operational opacity directory.

## 5. Running the program

### Check input files without loading opacities

```bash
python3 -m jupiter_rt --case . --preset co2 --check-inputs
```

This checks the `.fuel`, `.pta`, `.spe`, and instrumental-resolution inputs. It does **not** read `.opa` tables or run radiative transfer; it also does not independently verify the contents of `corrch4.txt`.

### Normal molecular abundance retrieval

```bash
python3 -m jupiter_rt \
  --case . \
  --opacity-dir /obs/tfouchet/data/opacite/jupiter/jwst \
  --preset abundance \
  --jacobian original \
  --output results_abundance
```

### Special-molecule retrieval: CO₂

```bash
python3 -m jupiter_rt \
  --case . \
  --opacity-dir /obs/tfouchet/data/opacite/jupiter/jwst \
  --preset co2 \
  --jacobian original \
  --output results_co2
```

For benzene, the presets `c6h6` and `c6h6_fortran` use different regularization settings. To select a different molecule, **change the input case files**, not merely the `--preset` argument.

### Forward calculation only

```bash
python3 -m jupiter_rt \
  --case . \
  --opacity-dir /obs/tfouchet/data/opacite/jupiter/jwst \
  --preset temperature \
  --forward-only \
  --output results_forward
```

### Long-running processes

```bash
nohup python3 -u -m jupiter_rt \
  --case . \
  --opacity-dir /obs/tfouchet/data/opacite/jupiter/jwst \
  --preset co2 \
  --output results_co2 \
  > log_co2.txt 2>&1 &
```

Monitor output with `tail -f log_co2.txt`; inspect processes with `ps -u "$USER" -o pid,etime,cmd | grep '[j]upiter_rt'`.

### Command-line arguments

| Argument | Purpose |
|---|---|
| `--case PATH` | Directory containing the input files; **required**. |
| `--opacity-dir PATH` | Directory containing `.opa` tables; defaults to `--case` when omitted. |
| `--preset NAME` | Built-in configuration; default `auto`. |
| `--jacobian METHOD` | `original` (default), `analytical`, or `finite_difference`. |
| `--output PATH` | Output directory, default `resultados_python`. Relative paths are relative to the current working directory. |
| `--check-inputs` | Check inputs without computing opacities. |
| `--forward-only` | Compute a synthetic spectrum without performing a retrieval. |
| `--options FILE.json` | Override retrieval options using a JSON file. |
| `--opacity-metadata FILE.json` | Specify alternate opacity-table grids. |
| `--endian little\|big` | Byte order of `.opa` data; default `little`. |
| `--max-samples N` | Spectral sample limit; default 160,000. |
| `--without-ch4-legacy-scale` | Disable the historical `7/6` scaling of the first gas. |

### Custom retrieval options

`--options` accepts a JSON object, for example:

```json
{
  "max_iterations": 30,
  "rtol_chi2": 0.001,
  "c_temperature": 0.75,
  "c_special": 100.0,
  "factor": 3.0,
  "stop_mode": "combined"
}
```

Only the specified settings override the selected preset. The `c_*` values control vertical correlation lengths, while `factor` enters the covariance scaling. **Retrievals using different regularization settings should not be compared as if every difference came from the forward model or Jacobians.**

## 6. Physical and numerical model

The implemented processing sequence is:

1. `io.py`: load pressure, temperature, gas profiles, and observations.
2. `geometry.py`: calculate altitude, latitude-dependent gravity, and viewing geometry.
3. `opacity.py`: read binary molecular opacities and interpolate them in temperature.
4. `cia_python.py` and `continuum.py`: add **H₂–H₂, H₂–He, and H₂–CH₄** collision-induced absorption (CIA). No external Fortran library is required during normal operation.
5. `forward.py`: accumulate optical depths and transmission, compute layer thermal emission, convolve with the instrument response, and calculate Jacobians when needed.
6. `inverse.py`: apply regularization, calculate parameter increments, and iteratively update temperature, VMR, and/or aerosols.

Monochromatic thermal radiance is constructed from the Planck function and differences in layer transmission. Each layer temperature is the average of the two adjacent level temperatures. The physical CIA formulation has not been replaced with a faster approximation that would alter the inherited numerical results.

Instrumental convolution is Gaussian, with an effective width controlled by the preset's `gaussian_width` parameter and the spectral resolution in `resolution.txt`.

### Opacity grids

Default coverage in `opacity.py`:

| Opacity-table type | Binary-block coverage |
|---|---|
| Normal gases | 500–2000 cm⁻¹ |
| HCN | 680–750 cm⁻¹ |
| CO₂ | 640–710 cm⁻¹ |
| H₂O | 500–2000 cm⁻¹ |
| C₆H₆ | 640–710 cm⁻¹ |

Each block spans 10 cm⁻¹, with 10,000 spectral samples at 0.001 cm⁻¹ spacing. These values describe **binary opacity-table coverage**, not necessarily a case's observed spectral interval. If the real opacity tables have different grids, supply appropriate information using `--opacity-metadata`; the code does not automatically detect their physical structure.

## 7. Retrieval and regularization

The retrieval parameter blocks are determined by `.fuel`: temperature is retrieved if any gas has flag 1 or 3; the VMR of gases with flag 2 or 3 is retrieved; aerosols are retrieved when `cloud` is nonzero.

The retrieval minimizes a regularized spectral mismatch using a Gaussian correlation in log-pressure:

\[
S_{ij}(c) = \exp\!\left[-\frac{1}{2}\left(\frac{\ln p_i - \ln p_j}{c}\right)^2\right].
\]

The code assigns `c_temperature`, `c_molecule`, `c_special`, or `c_aerosol` to the respective blocks. It calculates an adaptive `alpha` that includes the `factor` option, observational errors, and Jacobian sensitivity. It then solves a covariance-based system in spectral space and obtains parameter increments (`delta`), estimated uncertainties, and averaging kernels.

**Molecular updates are logarithmic**: `VMR_new = VMR_old * exp(delta)`. Accordingly, a retrieved gas must have a strictly positive starting VMR. Temperature is updated additively; aerosol optical depth is updated additively with a lower limit of zero.

Legacy vertical filters use hyperbolic-tangent masks centered near atmospheric level index 90. Their exact behavior depends on `filter_mode`:

- `legacy_abundance` and `legacy_co2`: apply a tropospheric filter to the first parameter block and a stratospheric filter to the second.
- `legacy_aerosol` and `legacy_c6h6`: apply a tropospheric filter to the first block.
- `legacy_hcn`: multiply the first block's increment by `scale_hcn_delta` (default 0.1).
- `legacy_h2o` and `none`: apply no additional mask in the current version's `apply_filters` function.

The currently implemented transition width is **15 levels**. This has been retained in v0.4.4 for continuity with the scientific experiments; it should not be assumed to reproduce every historical Fortran filter exactly.

### Stopping criteria

- `fortran`: tests the absolute relative χ² change using the **current** χ² in the denominator.
- `combined`: additionally requires `max|delta|` to decrease and fall below `tol_step` (default 0.3).
- The iteration limit is preset-dependent.

The iteration log reports the **fractional chi-squared improvement**:

\[
\frac{\Delta\chi^2}{\chi_{\mathrm{prev}}^2}
=\frac{\chi_{i-1}^2-\chi_i^2}{\chi_{i-1}^2}.
\]

The short output label is `Δchi2/chi2_prev`. Positive values mean the spectral fit has improved; negative values mean chi-squared increased. The first iteration reports `N/A`. Log fields are separated by `|`, for example:

```text
iter  3 | chi2_mean=0.838144 | Δchi2/chi2_prev=+8.582e-01 | max|delta|=0.00213
```

This display-only convention does not change the `rtol_chi2` stopping check, whose denominator is the **current** chi-squared. The JSON field `relative_chi2_change` uses the same positive-for-improvement convention.

The CLI also prints the **actual correlation length of each retrieved parameter** rather than a generic `c(mol)` value. For instance, `--preset co2` now displays `c(T)=0.75 | c(co2)=100`; this is only a reporting correction (the retrieval already used these values).

## 8. Available presets

Values defined by `cli.py` and `RetrievalOptions` in v0.4.4:

| Preset | `gaussian_width` | Relevant correlation lengths | Filter | Stopping mode |
|---|---:|---|---|---|
| `temperature` | 1.000 | `c_T=0.75` | `none` | `fortran` |
| `abundance` | 0.925 | `c_T=3`, `c_mol=3` | `legacy_abundance` | `fortran`, max. 20 |
| `abundance_split` | 0.925 | `c_T=0.75`, `c_mol=3` | `legacy_abundance` | `fortran`, max. 20 |
| `aerosol` | 1.000 | `c_aerosol=0.75` | `legacy_aerosol` | `combined` |
| `hcn` | 1.000 | `c_special=100` | `legacy_hcn` | `combined` |
| `co2` | 1.200 | `c_special=100` | `legacy_co2` | `combined` |
| `h2o` | 1.000 | `c_special=100` | `legacy_h2o` | `combined` |
| `c6h6` | 0.925 | `c_special=100` | `legacy_c6h6` | `combined` |
| `c6h6_fortran` | 0.925 | `c_T=2`, `c_special=2` | `legacy_c6h6` | `fortran`, max. 20 |
| `auto` | 0.925 | depends on retrieved parameter types | `auto` | `combined` |

Unless overridden explicitly, defaults are `c_T=0.75`, `c_mol=3`, `c_special=100`, `c_aerosol=0.75`, `factor=3`, `max_iterations=15`, and `rtol_chi2=0.01`. All listed presets use the overall `factor=3`.

**Important:** names such as `fortran` and `legacy` indicate compatibility goals, not bit-for-bit identity with every original Fortran program. The thermal Jacobian and some filters, for example, differ from certain historical variants.

## 9. Jacobians

The `--jacobian` option affects **molecular abundance Jacobians**, not the temperature or aerosol Jacobians.

- `original` (default): preserves the inherited expression using a spectrally subsampled grid (`nu[::5]`). It is the reference choice for continuity with existing results.
- `analytical`: evaluates an alternative molecular derivative based on transmission, using the full spectral grid.
- `finite_difference`: evaluates centered differences of **algebraic transmission expressions**, rather than independently running the full forward model twice with perturbed atmospheric states. It is an internal consistency check for the analytical formulation.

The thermal Jacobian distributes the layer Planck derivative across adjacent temperature levels. The aerosol Jacobian retains the project's original formulation, pending full independent validation.

The `analytical` and `finite_difference` molecular Jacobians were found to agree very closely for the studied abundance-retrieval case. That comparison **does not replace** a finite-difference test of the entire forward model. Numerical comparisons are documented in [VALIDATION.md](VALIDATION.md).

## 10. Reading and interpreting results

A retrieval creates these files under the directory specified by `--output`:

| File | Meaning |
|---|---|
| `jupiter_ch4.res` | Wavenumber, observed spectrum, synthetic spectrum, and measurement error (four columns). |
| `jupiter_ch4.inv` | Pressure, retrieved temperature, gas profiles, aerosols if active, and uncertainties; see its header. |
| `jupiter_ch4.sum` | Final iteration, mean χ², and DOFS for each retrieved parameter. |
| `jupiter_ch4.avg` | **1-based index of each averaging-kernel column's maximum** at each level; the full 361×361 matrix is not saved. |
| `diagnostics.json` | Final χ², iteration count, parameters, correlation lengths, Jacobian method, and iteration history. |

`--forward-only` writes spectra and profiles without the complete retrieval diagnostics.

The reported statistic is **mean χ²**, the mean squared residual normalized by observational error. It should not automatically be interpreted as a reduced χ² computed using an explicit residual degrees-of-freedom count.

DOFS are traces of the averaging-kernel matrices. The peak locations in `jupiter_ch4.avg` indicate where each column attains its maximum; they cannot, on their own, establish the kernel width, full shape, or number of vertically independent measurements.

**Scaling warning:** `io.py` applies the scale factors declared in `.fuel`, an additional historical `7/6` factor to the **first gas** unless `--without-ch4-legacy-scale` is set, and an additional **×50 only to C₆H₆**. For retrieved gases, `.inv` divides by the `.fuel` scale only, **not** by the historical `7/6` or `50` factors. Review the conventions before comparing `.inv` columns directly with raw `.pta` values.

## 11. Troubleshooting

| Symptom | Typical cause | What to check |
|---|---|---|
| `No module named jupiter_rt` | The package is not installed in the active interpreter. | Run `python3 -m pip install -e .` from the repository and test `python3 -c 'import jupiter_rt'`. |
| `FileNotFoundError: jupiter_ch4.fuel` | `--case` refers to the wrong directory. | Check `ls /path/to/case/jupiter_ch4.fuel`; `--case .` means the current working directory. |
| `c2h3` is retrieved after selecting `c6h6` | The input files describe another case. | Inspect `.fuel`, the `.spe` spectral range, and `Retrieved parameters` from `--check-inputs`. |
| Missing opacity file | The `.opa` file is absent or has a mismatched name. | Check `<opacity-dir>/jupiter_ch4_<gas>.opa`. |
| Unexpected opacity-table size | The actual binary layout differs from the configured one. | Verify the format and metadata; do not alter blocks without checking. |
| Extra `.pta` column | VMR and aerosol columns were confused. | The last column of a special-molecule case is its VMR; aerosols require no column. |
| Nonpositive molecular profile | Logarithmic retrieval requires positive initial VMR. | Check starting profiles for species with flags 2/3. |
| Large profile changes at similar χ² | Different regularization or weak vertical information content. | Compare presets, DOFS, uncertainties, and kernels, not just χ². |
| Missing `corrch4.txt` | The CIA calculation cannot be completed. | Place this file in the case directory and run an actual forward/retrieval case. |

## 12. Validation status and remaining work

- **CIA continuum:** individual components compared numerically against the original Fortran routines.
- **Temperature retrieval:** very close Python–Fortran agreement for the studied reference case.
- **Normal C₂H₂ abundance (`c2h3`) with thermal adjustment:** very close agreement with the reference retrieval.
- **C₆H₆:** the correct case runs and the modeled spectrum agrees very closely with Fortran; the vertical VMR differences are associated with different regularization settings. The current implementation is accepted for continued work, **without claiming identical retrieved profiles**.
- **CO₂:** compared successfully against Fortran for the shared special-molecule retrieval path. The 28-point spectra and stratospheric VMR agree well, although the Fortran and Python runs use different regularization and observation-error scaling.
- **HCN and H₂O:** not individually compared; a successful CO₂ test would not automatically verify their specific opacity tables, filters, or spectra.
- **Aerosols:** the inherited treatment is preserved, but a fully documented end-to-end comparison remains outstanding.
- **Jacobians:** `analytical` and `finite_difference` are internally consistent; full forward-model perturbation and comparison of complete averaging-kernel matrices are still pending.

Numerical results, experimental conditions, and interpretation caveats appear in [VALIDATION.md](VALIDATION.md). Future reference runs should record code version, preset, `.fuel`, spectral grid, regularization, Jacobian, iteration count, and Fortran baseline to avoid comparing mismatched cases.
