# Scientific Validation Status — Python versus Fortran

**Documented software:** Jupiter RT Python v0.4.4.  
**Status as of:** November 9, 2026.

## 1. Comparison principles

Reproducibility is assessed separately for:

1. **Forward model:** synthetic spectrum calculated from a prescribed atmosphere.
2. **Retrieval:** final atmospheric profiles and iterative convergence path.
3. **Information content:** degrees of freedom for signal (DOFS), estimated errors, and averaging kernels.
4. **Configuration:** identical case, observed spectrum, instrument response, `c_*`, `factor`, vertical filters, and stopping conditions.

A similar spectral fit does not, by itself, prove that an inverse solution is unique. Different regularization settings can yield different vertical profiles with almost identical χ².

The `.opa` files, observations, and large Fortran reference outputs are **not distributed with this source-code package**. The figures below summarize comparisons carried out during development; without the corresponding inputs, they are not automatically reproducible third-party tests.

## 2. Collision-induced absorption (CIA)

The Python continuum implementation includes H₂–H₂, H₂–He, and H₂–CH₄ without executing the original Fortran. Individual components were compared with outputs from the Fortran reference routines.

| Component | Largest observed relative difference |
|---|---:|
| H₂–H₂ | 5.644 × 10⁻⁷ |
| H₂–He | 3.322 × 10⁻⁷ |
| H₂–CH₄ | 3.411 × 10⁻⁷ |

The initial translation tests covered temperatures from 50 to 757.1 K and wavenumbers from 500 to 2000 cm⁻¹. An additional comparison of continuum optical depth at selected atmospheric levels found relative differences below 2.7 × 10⁻⁷.

**Scope:** These results establish numerical agreement for the tested conditions; they do not independently validate physical extrapolation beyond the applicability ranges of the original approximations.

## 3. Temperature retrieval

Historical reference case with **361 levels and 373 spectral samples**:

| Metric | Fortran | Python |
|---|---:|---:|
| Final mean χ² | ≈ 3150.2481 | ≈ 3149.9217 |
| Temperature RMS difference | — | ≈ 0.028 K |

The spectral and thermal discrepancies between implementations are small for the studied case. Python uses a thermal Jacobian distributed between neighboring temperature levels, which is not a literal reproduction of every detail in the original formulation.

## 4. Normal molecular abundance retrieval: C₂H₂ (`c2h3`)

Validated case with tropospheric temperature adjustment and abundance retrieval, **361 levels and 133 points** spanning 685.13–719.93 cm⁻¹:

| Metric | Fortran | Python |
|---|---:|---:|
| Mean χ² | ≈ 0.229879 | ≈ 0.229197 |
| Final iteration | 9 | 9 |
| Temperature RMS difference | — | ≈ 0.0173 K |
| Relative abundance RMS difference | — | ≈ 0.128% |

This comparison used the historical abundance settings **`c_temperature=3`, `c_molecule=3`, `factor=3`**, with `--jacobian original`. The gas and opacity identifier remains `c2h3` even when the molecule is scientifically discussed as C₂H₂.

**Conclusion:** Excellent numerical agreement for the tested case; this does not independently validate the inherited molecular Jacobian formula.

## 5. Molecular Jacobian comparison

The `original`, `analytical`, and `finite_difference` methods were compared using the normal abundance case:

| Method | Final mean χ² | Final iteration | Temperature DOFS | Gas DOFS |
|---|---:|---:|---:|---:|
| `original` | ≈ 0.229197 | 9 | ≈ 1.01186 | ≈ 1.41289 |
| `analytical` | ≈ 0.231308 | 8 | ≈ 1.00608 | ≈ 1.41258 |
| `finite_difference` | ≈ 0.231308 | 8 | ≈ 1.00608 | ≈ 1.41258 |

Part of the difference in final χ² arises from **different stopping iterations**. At iteration 8, the original method had χ² ≈ 0.231327, essentially the same as the analytical method.

`analytical` and `finite_difference` agree to very high numerical precision. However, the latter uses centered differences of **an algebraic transmission expression** and shares code with the analytical formulation. **An independent perturbation of the complete forward model has not yet been performed** to test these methods against the full physical derivative. The original Jacobian also uses a spectrally subsampled grid, whereas the alternative methods use the full grid; this distinction can contribute to the observed differences.

The vertical distributions of peak positions in `jupiter_ch4.avg` are not identical, particularly for the gas. Similar DOFS do not establish identical averaging-kernel shapes. That would require exporting and comparing the complete averaging-kernel matrices.

**Current decision:** Retain `original` as the default molecular Jacobian.

## 6. C₆H₆: evaluated special-molecule case

The correctly matched Fortran `00-01` case was compared with its corresponding Python run:

- **361 levels and 54 points** spanning 671.11–676.94 cm⁻¹.
- Retrieved parameters: temperature and C₆H₆.
- Both runs ended at **iteration 8**.

| Metric | Fortran | Python |
|---|---:|---:|
| Final mean χ² | ≈ 0.0803665 | ≈ 0.080999 |
| Temperature DOFS | ≈ 0.954 | ≈ 0.95416 |
| C₆H₆ DOFS | ≈ 0.953 | ≈ 0.95273 |

The relative difference in mean χ² is approximately **0.79%**. The maximum difference between synthetic spectra was about **1.32 × 10⁻¹¹**, corresponding to **0.0085 times** the relevant observational uncertainty. The temperature RMS difference was about **0.216 K**, with a maximum difference of approximately **0.814 K**.

**Caveat:** The retrieved C₆H₆ VMR profiles exhibited noticeable vertical differences (factors of two or three in some regions), despite excellent spectral agreement. According to the comparison of the configurations, the runs used **different regularization values**, so these discrepancies should not be attributed to an implementation error. The scientific decision was to **keep v0.4.4 unchanged**.

This result supports operation of the C₆H₆ forward model and special-molecule retrieval path. It **should not be presented as an exact profile-by-profile Fortran reproduction under identical mathematical settings**.

### Special scaling

The code reads the C₆H₆ VMR from the **last `.pta` column**, applies the general scaling conventions, and **multiplies it by an additional factor of 50 only for `c6h6`**, reproducing an identified historical feature of the reference Fortran program. The ×50 factor does not apply to HCN, CO₂, or H₂O.

## 7. CO₂: evaluated special-molecule case

CO₂ was compared with a historical Fortran retrieval using **361 atmospheric levels and 28 spectral samples** between 665.93 and 668.93 cm⁻¹. Both calculations retrieved temperature and CO₂. They used different settings; in particular, historical Fortran applied `profil(:,7) *= 50` to the CO₂ prior and multiplied its observational uncertainties by 75, while Python does **not** apply either CO₂-specific multiplier. The Fortran output uses a total chi-squared and Python reports the mean per spectral point.

| Comparison metric | Result |
|---|---:|
| Fortran total χ² using Fortran errors | ≈ 0.046201 |
| Python total χ² recomputed with Fortran errors | ≈ 0.047197 |
| Temperature RMS difference | ≈ 0.185 K |
| Temperature maximum absolute difference | ≈ 0.50 K |
| CO₂ mean absolute VMR difference, 10⁻⁸–10⁻² bar | ≈ 2.20% |
| CO₂ maximum absolute VMR difference, 10⁻⁸–10⁻² bar | ≈ 4.18% |
| Maximum synthetic-spectrum difference | ≈ 3.97 × 10⁻¹¹ |
| Maximum spectral difference in Fortran uncertainty units | ≈ 0.00263 σ |

The Fortran and Python temperature and CO₂ DOFS were also similar (roughly 0.92 and 0.91, respectively). This is **good agreement for the shared special-molecule calculation and the stratospheric CO₂ profile**, not a claim of exact equivalence of every posterior quantity. The temperature regularization, special-molecule regularization, vertical filter widths, and observational errors differ between the historical runs. Some Fortran uncertainty outputs contained non-finite values, so the posterior error estimates are not yet fully cross-validated.

## 8. Aerosols and remaining special molecules

| Component | Current status | Remaining validation |
|---|---|---|
| Aerosols | Legacy retrieval treatment preserved; initial aerosol profile is zero. | End-to-end Fortran comparison, including masks, spectra, retrieved opacity, and averaging kernels. |
| HCN | Shared special-molecule infrastructure tested indirectly through CO₂; HCN not tested individually. | Validate HCN opacity data, spectral interval, and case-specific settings if needed. |
| H₂O | Shared special-molecule infrastructure tested indirectly through CO₂; H₂O not tested individually. | Validate H₂O opacity data, spectral interval, and case-specific settings if needed. |

The **main remaining retrieval family to test is aerosols**. HCN and H₂O are intentionally deferred individual checks, not separately validated cases. The original Fortran sources for both molecules have also been inspected: they apply `profil(:,7) *= 1` to the special-molecule profile and do not multiply observational errors. Their only prior scaling consists of the shared `.fuel` factors and the first-gas `1.166667` correction. This supports the current Python input-scaling convention for HCN and H₂O but **does not replace full numerical tests** of their opacity files and retrievals.

A complete forward-model finite-difference Jacobian test and exports of full averaging-kernel matrices are possible future refinements, not changes introduced in v0.4.4.

## 9. Limits of the validation claims

These results represent specific test cases, not a general reproducibility certificate for every latitude, species, instrumental resolution, or thermodynamic range. For scientific publications, retain the code version and retrieval parameters, and rerun the relevant comparisons using the actual inputs and numerical precision of the intended analysis.
