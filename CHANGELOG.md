# Changelog

## 1.2.1 — 2026-11-13

- Corrected forbidden-transition decoding in the beta-minus and beta-plus engines. Unique forbidden shapes use the actual spin change; non-unique and unresolved transitions use an explicit allowed-shape approximation.
- Replaced obsolete signed forbiddenness codes with shared spin/parity transition metadata and simplified branch editing.
- Batched endpoint Monte Carlo calculations while preserving independent samples and the previous random draw order.
- Fixed compiled BSG array broadcasting and beta-plus sign selection, and improved CMake compatibility.
- Regenerated the default ENSDF and XUNDL neutrino spectra with the corrected transition handling.
- Included the shared transition module and BSG subpackage in package discovery, with regression tests for transitions and endpoint uncertainty.

## 1.2.0 — 2026-10-22

- Added ENSDF decay-mode normalization (`BR`/`NB`) and isotope-level `beta_frac` support.
- Added conservative missing-branch completion, with optional XUNDL supplementation matched by daughter excitation energy.
- Added `BetaEngine(process_missing=..., use_xundl=...)` configuration; summation uses the data selected when the beta engine is built.
- Regenerated the packaged ENSDF beta-minus, EC/beta-plus, and XUNDL databases.
- Added tentative/provenance metadata to parsed ENSDF isotope and branch records.
- Verified packaged ENDF fission-product yields against the source ENDF/B-VIII.0 data and `FPYEngine` loading.
