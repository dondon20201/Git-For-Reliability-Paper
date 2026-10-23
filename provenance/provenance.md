# Provenance and reproducibility

This is from the development repository, from the file `analysis/e11_manifest.R`. The dates below are from that repository's history. Because these commits are not in this repository's Git history, the order of these commits cannot be verified here.

## Primary data: CPC18 estimation raw data

- Source: Zenodo record 845873 (CC-BY), fetched and verified by `R/fetch_cpc18.R`.
- Raw CSV md5: `38a26407a36c690be969f9d9d3fe3e5a`
- Data-dictionary PDF md5: `3980a6da798b4883167ebfcc7a2c52c5`
- The raw data is NOT redistributed in this repository (gitignored). Reproducibility comes from the verified fetch script plus the renv lockfile. On an md5 mismatch the fetch stops and nothing is analysed.

## Secondary data: Basel-Berlin Risk Study (Frey et al.)

- Source: OSF project rce7g (Basel-Berlin Risk Study, Frey et al. 2017, Science Advances); local copy at data/raw/Frey Data/ (gitignored, never committed)
- Checksums for the 15 Frey input files are listed in `data/README.md`. Step 27 writes `output/elevation/results/p3_frey/provenance_manifest.csv` at runtime and also records `Codebook.pdf` for documentation.
- Not redistributed here; the local copy is gitignored.
- Used for external replication only (P3). It cannot and does not support any CPC18 per-subject claim.

## Pre-specification timeline

These dates are from the history of the development repository. They are preserved here for archival purposes, and are not dates defined by commits within this repository:

| Specification | Spec committed | Results committed | Order |
|---|---|---|---|
| `elevation/E3_spec.md` | 2026-07-17 | 2026-07-18 | spec first |
| `elevation/E3b_spec.md` | 2026-07-18 | 2026-07-18 | spec first |
| `elevation/E4_spec.md` | 2026-07-18 | 2026-07-18 | spec first |
| `elevation/E8_spec.md` | 2026-07-19 | 2026-07-19 | spec first |
| `elevation/P2_spec.md` | 2026-07-19 | 2026-07-20 | spec first |
| `elevation/P3_spec.md` | 2026-07-19 | 2026-07-19 | spec first |
| `elevation/SENS_spec.md` | 2026-07-21 | 2026-07-21 | spec first |
| `elevation/DISPB_spec.md` | 2026-07-22 | 2026-07-22 | spec first |
| `elevation/DISPC_spec.md` | 2026-07-22 | 2026-07-22 | spec first |

Pre-declared thresholds were fixed before results existed; the pre-declared failures are reported as results (see `elevation/results/hypothesis_ledger.csv`).

## Environment

- R 4.6.1, dependencies pinned by `renv.lock` (`renv::restore()` reproduces the library).
- In the development repository, the manuscript was generated using Quarto to pdf (typst) and docx. The numerical values were rendered using `num()` from the file `results/step5/manuscript_numbers.csv`. A render check was implemented in the file `R/verify_manuscript_numbers.R`. The manuscript files and render check are not included in this repository. The preprint was then modified after this step.

## Analysis manifest

`../manifest/analysis_manifest.csv` maps the paper's figures, tables, and result groups to analysis steps and results files. It explicitly marks a file absent from `expected/` rather than implying that it is committed.

