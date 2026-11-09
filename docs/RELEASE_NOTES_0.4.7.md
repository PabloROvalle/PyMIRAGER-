# PyMIRAGER v0.4.7 — release notes

**Release date:** 2026-11-10  
**Maturity:** 0.x scientific research release  
**Target:** Jupiter, thermal-infrared emission in local thermodynamic equilibrium (LTE)

## Scope

This release implements temperature and selected molecular retrievals for the existing Jupiter opacity-table/case conventions, and **grey absorbing aerosols**. Other planets, unspecified molecules, aerosol scattering, wavelength-dependent aerosol extinction, and non-LTE physics are not implemented. They are potential future projects, not existing features.

## The one aerosol-physics configuration change

`--preset aerosol_cutoff` now sets:

```json
{"aerosol_zero_below_index": 85, "aerosol_zero_from_index": 120}
```

Consequently `taucloud[0:85]` and `taucloud[120:361]` are zero, and **indices 85–119** can contain retrieved optical depth. Both aerosol Jacobian sensitivity and retrieval updates are zeroed outside this band. The separate cloud-mask file continues to apply. There is no scattering and no new pressure weighting; `taucloud` is a layer optical depth.

The previously assessed run used **80–119**, not 85–119. Its case-specific results (mean χ² approximately 6.696821; integrated tau approximately 0.004221) **must not be quoted as new 85–119 results**. A fresh full retrieval is required to validate the new boundary scientifically.

## Additional changes

- Unified package and distribution version `0.4.7`.
- Added `python3 -m jupiter_rt --version` and extra provenance in `diagnostics.json`.
- Reworked documentation, added version history and offline unit tests.
- Preserved forward physics, numerical inversion formulae, all other preset settings, and data format compatibility.

## Local verification

```bash
python3 -m pip install -e .
python3 -m unittest discover -s tests -v
python3 -m jupiter_rt --version
cd ~/MIRAGER
python3 -m jupiter_rt --case . --preset aerosol_cutoff --check-inputs
```

For the actual science run, provide the external `.opa` tables and mask and check that the diagnostics report `aerosol_zero_below_index: 85` and `aerosol_zero_from_index: 120`.

## Publishing

After reviewing files and verifying the repository's working tree:

```bash
git add README.md CHANGELOG.md pyproject.toml jupiter_rt/cli.py jupiter_rt/__init__.py docs tests
git commit -m "Release PyMIRAGER v0.4.7 (Jupiter LTE scope, aerosol 85-119)"
git tag -a v0.4.7 -m "PyMIRAGER v0.4.7"
git push origin HEAD
git push origin v0.4.7
```

These are instructions; preparing the source bundle **does not by itself publish a GitHub release or tag**.
