# Changelog

Changes are grouped by release. All versions remain in the **0.x research series**; a 1.0 general-purpose planetary interface has not been defined.

## [0.4.7] — 2026-11-10

### Changed
- `aerosol_cutoff`: move the lower fixed-zero boundary from **index 80 to 85**. The allowed zero-based vertical indices are **85–119**; the upper boundary remains **120**.
- Set the distribution version in `pyproject.toml` and `jupiter_rt.__version__` consistently to **0.4.7**.
- Reformat preset definitions for maintainability (all other preset values unchanged).
- Add `--version` to the CLI and include `software_version` and `preset` in retrieval `diagnostics.json`.
- Update the README and reference documentation to state the actual **Jupiter-only, LTE, 5–20 µm, supported-opacity, grey-absorption aerosol** scope and the future roadmap.
- Record the already completed grey-aerosol comparison; distinguish the previously tested **80–119** band from the new **85–119** configuration, which requires a fresh end-to-end scientific run.

### Added
- Offline unit/regression tests for the preset, hard-zero aerosol constraints, version consistency, and invalid bounds.
- Release notes and reproducible local release instructions.

### Unchanged
- Thermal radiative transfer, CIA, gas opacities, regularization formulae, molecule Jacobians, and other retrieval presets.
- HCN and H₂O remain individually unvalidated.

## [0.4.6] — 2026-11-09
- Add a second aerosol cutoff, allowing levels 80–119 in `aerosol_cutoff`.
- Compare that configuration with Fortran for the documented reference case.

## [0.4.5] — 2026-11-09
- Add a selectable upper-only aerosol cutoff at index 120.

## [0.4.4] — 2026-11-09
- Modular Jupiter forward/retrieval implementation and scientific comparisons documented in `docs/VALIDATION.md`.
