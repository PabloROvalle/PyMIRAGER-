# PyMIRAGER
### Python-based Multi-target InfraRed Atmospheric Giant-planet Emission Retrieval

🚧 **Under active development**

## Overview

**PyMIRAGER** is a Python-based atmospheric retrieval framework designed to investigate the thermal structure, chemical composition, and aerosol properties of the four giant planets of the Solar System: Jupiter, Saturn, Uranus, and Neptune.

The framework is specifically conceived for the analysis of mid-infrared thermal emission spectra acquired with the **Mid-Infrared Instrument (MIRI)** aboard the **James Webb Space Telescope (JWST)**, particularly its **Medium Resolution Spectrometer (MRS)**.

PyMIRAGER solves the radiative transfer equation under the assumption of **Local Thermodynamic Equilibrium (LTE)**, enabling the simulation of emergent thermal emission spectra from giant-planet atmospheres.

The retrieval framework is designed to constrain:

- **Atmospheric temperature:** Vertical temperature profiles as a function of pressure.
- **Chemical composition:** Molecular abundances and their vertical distributions.
- **Aerosol properties:** The contribution of atmospheric aerosols to mid-infrared thermal emission spectra.

Although PyMIRAGER is intended to become a multi-target retrieval framework for all four giant planets, **the current implementation is configured for Jupiter**. Future developments will extend its applicability to Saturn, Uranus, and Neptune.

## Molecular Opacity and Continuum Absorption

PyMIRAGER incorporates molecular absorption and collision-induced absorption (CIA) processes relevant to giant-planet atmospheres.

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

The radiative transfer calculations require precomputed, species-specific opacity files (`.opa`), which can be provided upon request.

## Validation and Current Development Status

PyMIRAGER is currently undergoing development, testing, and validation of its atmospheric retrieval capabilities using mid-infrared observations of Jupiter.

**Successfully tested:**

- **Temperature retrievals:** The retrieval of vertical atmospheric temperature profiles has been successfully tested, yielding very promising results.
- **Acetylene (C₂H₂) retrievals:** The retrieval of acetylene abundances has also been successfully tested, demonstrating very good performance.

**Ongoing development and validation:**

- Testing and validation of retrievals for the remaining molecular species.
- Implementation and validation of aerosol-related retrieval capabilities.
- Further assessment of retrieval performance under different atmospheric conditions.
- Extension of the framework to Saturn, Uranus, and Neptune.

The retrieval methodologies and ongoing developments build upon the approaches presented in **Rodríguez-Ovalle et al. (2024a, 2024b, 2025)**.

While temperature and C₂H₂ retrievals have already demonstrated promising performance, the remaining molecular and aerosol retrieval capabilities are still under evaluation.

## Opacity Files and Availability

PyMIRAGER requires precomputed opacity files (`.opa`) containing the spectroscopic information needed for the radiative transfer calculations.

These files are not included in the standard distribution but can be made available upon request.

## References

1. Rodríguez-Ovalle et al. (2024a). *Temperature and Composition Disturbances in the Southern Auroral Region of Jupiter Revealed by JWST/MIRI*. Journal of Geophysical Research: Planets, 129, e2024JE008299.

2. Rodríguez-Ovalle et al. (2024b). *Stratospheric aerosols and C₆H₆ in Jupiter's south polar region from JWST/MIRI observations*. Astronomy & Astrophysics, 691, A51.

3. Rodríguez-Ovalle et al. (2025). *JWST observations of exogenic species on Jupiter: HCN, H₂O, and CO₂*. Astronomy & Astrophysics, 696, A173.
