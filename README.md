# PyMIRAGER
### Python-based Multi-target InfraRed Atmospheric Giant-planet Emission Retrieval

**PyMIRAGER** is a PYTHON BASED mid-infrared atmospheric retrieval framework designed to investigate the thermal structure, chemical composition, and aerosol properties of the four Solar System giant planets: **Jupiter, Saturn, Uranus, and Neptune**.

The code is primarily developed for the analysis of thermal emission spectra acquired with the **Mid-Infrared Instrument (MIRI), Medium Resolution Spectrometer (MRS), aboard the James Webb Space Telescope (JWST)**.

> 🚧 **Development status: Under active development and validation.**

## Overview

MIRAGER solves the radiative transfer equation under the assumption of **Local Thermodynamic Equilibrium (LTE)** to simulate the emergent thermal emission from giant-planet atmospheres in the mid-infrared spectral region.

By comparing synthetic spectra with observations, the framework aims to retrieve key atmospheric properties, including:

- **Temperature profiles:** Vertical distributions of atmospheric temperature as a function of pressure.
- **Molecular abundances:** Atmospheric composition and the vertical distribution of infrared-active species.
- **Aerosol properties:** Contributions of atmospheric aerosols to the observed infrared emission spectra.

The multi-target architecture is intended to provide a consistent framework for atmospheric retrievals across the four giant planets, while accounting for their individual atmospheric characteristics.

## Molecular Opacity and Continuum Absorption

MIRAGER incorporates molecular absorption and collision-induced absorption (CIA) processes relevant to giant-planet atmospheres.

**Continuum opacity:**
- H₂–H₂ collision-induced absorption
- H₂–He collision-induced absorption
- H₂–CH₄ collision-induced absorption

**Molecular species:**
- CH₄ — Methane
- CH₃D — Monodeuterated methane
- NH₃ — Ammonia
- PH₃ — Phosphine
- C₂H₂ — Acetylene
- C₂H₆ — Ethane
- C₆H₆ — Benzene
- HCN — Hydrogen cyanide
- CO₂ — Carbon dioxide
- H₂O — Water

The radiative transfer calculations require species-specific opacity files (`.opa`), which are not included in the standard distribution and can be provided upon request.

## Validation and Current Development Status

MIRAGER is currently undergoing testing and validation of its atmospheric retrieval capabilities.

**Successfully tested:**
- **Temperature retrievals:** The retrieval of atmospheric temperature profiles has been tested, yielding very promising results.
- **Acetylene (C₂H₂) retrievals:** The retrieval of acetylene abundances has also been successfully tested, showing very good performance.

**Ongoing validation:**
- Retrievals of the remaining molecular species supported by the framework.
- Implementation and testing of aerosol-related retrieval capabilities.
- Further evaluation of retrieval performance across different atmospheric conditions and planetary targets.

The development and validation of these capabilities build upon the radiative transfer and atmospheric retrieval methodologies presented in **Rodríguez-Ovalle et al. (2024a, 2024b, 2025)**.

Although the initial temperature and C₂H₂ retrieval tests have produced encouraging results, the remaining molecular and aerosol retrieval capabilities are still being evaluated.

## Opacity Files and Availability

MIRAGER requires precomputed opacity files (`.opa`) containing the spectroscopic information necessary for the radiative transfer calculations.

These files are available upon request. Users interested in obtaining the required opacity data or additional information about the framework are encouraged to contact the developer.

## References

The development of MIRAGER builds upon previous studies of giant-planet atmospheric retrievals using JWST/MIRI observations:

1. Rodríguez-Ovalle et al. (2024a). *Temperature and Composition Disturbances in the Southern Auroral Region of Jupiter Revealed by JWST/MIRI*. Journal of Geophysical Research: Planets, 129, e2024JE008299.

2. Rodríguez-Ovalle et al. (2024b). *Stratospheric aerosols and C₆H₆ in Jupiter's south polar region from JWST/MIRI observations*. Astronomy & Astrophysics, 691, A51.

3. Rodríguez-Ovalle et al. (2025). *JWST observations of exogenic species on Jupiter: HCN, H₂O, and CO₂*. Astronomy & Astrophysics, 696, A173.
