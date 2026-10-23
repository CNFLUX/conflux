# Changelog

## 1.2.0 — 2026-10-22

- Added ENSDF decay-mode normalization (`BR`/`NB`) and isotope-level `beta_frac` support.
- Added conservative missing-branch completion, with optional XUNDL supplementation matched by daughter excitation energy.
- Added `BetaEngine(process_missing=..., use_xundl=...)` configuration; summation uses the data selected when the beta engine is built.
- Regenerated the packaged ENSDF beta-minus, EC/beta-plus, and XUNDL databases.
- Added tentative/provenance metadata to parsed ENSDF isotope and branch records.
- Verified packaged ENDF fission-product yields against the source ENDF/B-VIII.0 data and `FPYEngine` loading.
