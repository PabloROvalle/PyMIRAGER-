# Input Data Formats, Scaling Factors, and Output Files

**Version:** Jupiter RT Python 0.4.4. This document describes the format **implemented by the current input reader**, not a generic specification for every historical Fortran input file.

## 1. Case directory layout

```text
case_directory/
  jupiter_ch4.fuel
  jupiter_ch4.pta
  jupiter_ch4.spe
  resolution.txt
  corrch4.txt
  nh3_clouds.txt     # required only for cloud 1
  ch4_hazes.txt      # required only for cloud 2
```

Binary `.opa` files are usually stored separately and located through `--opacity-dir`.

## 2. `jupiter_ch4.fuel`

The reader removes blank lines and text following `!` or `#` comments. The remaining records must appear in this order:

1. Latitude, in degrees.
2. Cosine of the incidence angle.
3. Cosine of the emergence angle.
4. `cloud N`, where `N` is 0, 1, or 2.
5. One line per gas with three fields: **gas name, scale factor, retrieval flag**.

Illustrative example for a special-molecule case:

```text
1.003336
0.994964
0.984629
cloud 0
ch5 1.0 1
ch3e 1624.0 0
nh3 1.0 0
ph3 1.0 0
c2h3 1.0 0
c2h6 1.0 0
c6h6 1.0 2
```

The numbers only illustrate the ordering; use scientifically appropriate values from the observation being analyzed. Gas identifiers correspond to opacity filenames and are not necessarily standard chemical formulas.

### Retrieval flags

| Flag | Meaning |
|---:|---|
| 0 | Fixed gas: its profile and opacity are used, but its VMR is not retrieved. |
| 1 | Enables **temperature retrieval** because at least one gas has flag 1; the gas VMR itself remains fixed. |
| 2 | Retrieves the gas **VMR**. |
| 3 | Enables both temperature retrieval and retrieval of this gas's **VMR**. |

Temperature is represented by a single parameter block even if several gases have flags 1 or 3. Choosing a preset alone does not activate a molecule: the relevant `.fuel` flag must also be set.

If a special molecule is included, it **must be the last gas** in the list.

### Aerosol modes

- `cloud 0`: no aerosol retrieval.
- `cloud 1`: retrieve aerosols using `nh3_clouds.txt`.
- `cloud 2`: retrieve aerosols using `ch4_hazes.txt`.

In every mode, `taucloud` is initially zero. The cloud files supply the **vertical mask** (361 rows containing pressure and mask), not an aerosol-opacity column in `.pta`.

## 3. `jupiter_ch4.pta`

This is an ASCII numeric table containing **361 rows**:

```text
pressure_bar  temperature_K  VMR_gas_1  VMR_gas_2  ...  VMR_gas_N
```

Column 1 contains decreasing pressure in bar; column 2 contains temperature in K. The subsequent columns follow the exact gas order in `.fuel`.

Example **column order**, not actual data values:

```text
P | T | ch5 | ch3e | nh3 | ph3 | c2h3 | c2h6 | c6h6
```

**Permanent rule:** there is NO aerosol column. For HCN, CO₂, H₂O, and C₆H₆, **the final column contains that molecule's VMR**.

For backward compatibility, v0.4.4 accepts one extra trailing column consisting entirely of zeros **only for cases without a special molecule**. That extra column is ignored; it is never interpreted as aerosol optical depth.

### Input scaling rules

The reader applies the following operations, in order:

1. Load all gas VMR profiles from `.pta`.
2. Multiply the **first molecular profile** by `7/6` for historical compatibility (disable this with `--without-ch4-legacy-scale`).
3. Multiply every gas profile by its numeric scale factor from `.fuel`.
4. If a gas is named exactly `c6h6`, multiply its profile **once more by 50**.

**HCN, CO₂, and H₂O do not receive the extra ×50 factor in the current Python code.** They do not receive an automatic ×5 factor simply because they are special molecules.

**Historical Fortran comparison:** the HCN and H₂O variants contain `profil(:,7) = profil(:,7)*1`, so neither adds a special-molecule prior multiplier beyond the generic `.fuel` and first-gas (`7/6`) scaling. Both variants also retain the spectral errors read from `.spe`, without a special multiplier. In contrast, the historical CO₂ variant multiplies the special-molecule prior by **50** and the spectral errors by **75**. These CO₂-specific Fortran operations are **not enabled in the current Python `co2` preset**: preserve this distinction when comparing old output files. The historical C₆H₆ and current Python `c6h6` implementations both apply the extra ×50 prior multiplier.

The `7/6` scaling applies to the first gas listed, even when its identifier is not literally `ch4`. Preserve the gas ordering expected by the original reference implementation.

### Interpreting VMR values in `.inv`

When writing retrieved gas profiles, the output routine divides by the scale factor specified in `.fuel`, but **does not undo** the additional historical `7/6` or `50` multipliers. Review the scaling convention before comparing a VMR column in `.inv` directly with raw values from `.pta`.

## 4. `jupiter_ch4.spe`

An ASCII file containing the number of observed points, `N`, in the first record, followed by **N rows** of:

```text
wavenumber_cm-1   radiance_observed   sigma_radiance
```

The row count must agree with `N`; uncertainties must be strictly positive. Radiance units are inherited from the original Fortran cases. **The implementation does not automatically convert arbitrary radiance units.**

## 5. `resolution.txt`

The first line specifies the number of tabulated rows. Each subsequent row contains a **wavenumber** and an instrumental **FWHM**, expressed in the model's wavenumber units.

The reader sorts the tabulated wavenumbers, linearly interpolates FWHM to each `.spe` sample, and rejects observed spectra outside the available resolution range. It does not extrapolate.

Instrumental convolution also depends on the preset's `gaussian_width`; this value must be considered when comparing against a Fortran run.

## 6. `corrch4.txt`

Contains the historical coefficient table for the H₂–CH₄ collision-induced absorption continuum, read by `cia_python.py`. The model also calculates H₂–H₂ and H₂–He CIA. Keep the coefficient file associated with the original case; the reader expects its historical layout. `--check-inputs` does not verify these coefficients, because CIA is evaluated during the forward calculation.

## 7. `.opa` opacity tables

One file is expected for each species in `.fuel`:

```text
<opacity-dir>/jupiter_ch4_<gas_name>.opa
```

These are binary `float32` tables with 361 vertical levels, five tabulated temperatures per level, and 10 cm⁻¹ spectral blocks. The storage order and encoding are part of the historical compatibility requirements. Do not convert the files to text or alter their layout without updating `opacity.py` and the corresponding grid definition.

The `--endian` option selects the byte order (`little` by default). The `--opacity-metadata` option accepts a JSON object keyed by gas name. Its values may define `OpacityGrid` fields: `origin`, `blocks`, `nfreq`, `dnu`, `ntemp`, and `nlevel`.

The complete opacity tables can be large and are not included in the public source-code repository. Opacity locations are configurable paths, **not hard-coded personal directories**.

## 8. Output files

### `jupiter_ch4.res`

Four columns: `nu_cm^-1`, `observed`, `synthetic`, `error`.

### `jupiter_ch4.inv`

- `p_bar`, `T_K`.
- One column per retrieved gas (flag 2 or 3), in `.fuel` order, expressed as `gas_profile_over_scale`.
- `tau_cloud` if aerosol retrieval is active.
- One estimated-uncertainty column, `sigma_<parameter>`, for each retrieved parameter.

The file includes a **commented header** identifying the columns used in that specific run. Analysis scripts should read this header rather than assume fixed positions.

### `jupiter_ch4.sum`

A single row containing the final iteration, mean χ², and `dof_T`, `dof_<gas>`, or `dof_aerosol` as appropriate for the retrieved parameters.

### `jupiter_ch4.avg`

One row per atmospheric level (1–361) and one column per retrieval parameter, containing the **1-based index at which the corresponding averaging-kernel column has its maximum**. This is a peak-location summary, **not the complete averaging-kernel matrix**.

### `diagnostics.json`

Contains `chi2_mean`, `iterations`, `parameters`, `correlation_length`, `filters`, `molecular_jacobian`, and iteration `history`. The field `relative_chi2_change` stores `(previous_chi2 - current_chi2) / previous_chi2`: positive values indicate an improvement (chi-squared decreased); negative values indicate deterioration. Its first value is `null` (reported as `N/A` on screen). This field changed sign in v0.4.4; the retrieval stopping criterion did not change.

## 9. Checklist for Python–Fortran comparisons

For **both runs**, verify: input case, number of spectral samples, spectral interval, gases and flags, a priori profile, scale factors, preset, `c_*` and `factor`, vertical filters, Jacobians, stopping criterion, and software version. Similar χ² values **do not imply** identical retrieved vertical profiles when the regularization settings differ.
