# Analysis materials

The analysis code and reference results for "How Many Problems Are Enough, and Which? A Reliability Ceiling for Individual-Difference Measures of Risky Choice" are in this repository. These folders hold the related documentation and checks:

| Folder | Contents |
|---|---|
| `specifications/` | Analysis specifications, errata, and the project's integrity rules |
| `provenance/` | A reported timeline of the specifications and results in the development history |
| `manifest/` | A map from paper figures, tables, and main results to steps and results files |
| `figures/` | The script for drawing Figures 1 to 3 and the three PDF figures |
| `number_check/` | A checker that looks for paper numbers in the results files |

## Specifications and analysis steps

| Specification | What it specifies | Analysis step |
|---|---|---|
| `E1_contract.md` | Targets, participant split, H1, and H2 | 20 |
| `E2_spec.md` | Problem-information dictionary | 21 |
| `E3_spec.md` | Equal-length subset benchmark (H1) | 22, 23 |
| `E3b_spec.md` | High-precision subset reliability (exploratory) | 24 |
| `E4_spec.md` | Mechanism model (H2) | 25 |
| `E4b_spec.md` | Exploratory re-test of the mechanism on E3b outcomes | No implementation or results file is included here |
| `E8_spec.md` | Hierarchical re-estimation | 26 |
| `P3_spec.md` | Basel-Berlin Risk Study re-analysis | 27 |
| `P2_spec.md` and `P2_spec_erratum.md` | Joint item-information model | 28 |
| `SENS_spec.md` and `SENS_spec_erratum.md` | Mechanism sensitivity and attenuation bound | 29 |
| `DISPB_spec.md` | Partition stability, selection optimism, leave-one-feature-out | 30 to 33 |
| `DISPC_spec.md` | Local dependence and response-style diagnostics | 34, 35 |

Steps 36 and 37 assemble summary values and the hypothesis table. See `specifications/RELEASE_NOTES.md` for changes made to the archival copies of the specifications.

## Provenance

The timeline in `provenance/provenance.md` reports dates from an earlier development history. The specifications were added to this repository later. Its Git history does not by itself establish the reported order of specifications and results. The specifications also retain references to the earlier file layout. In particular, manuscript source and render checks mentioned in the integrity rules are not included here.

## Figures

From the repository root, with Python 3, pandas, numpy, and matplotlib installed:

```
python figures/make_figures.py expected figures
```

The script reads the committed results and writes `fig1_reliability.pdf`, `fig2_projection.pdf`, and `fig3_prior_studies.pdf`. Use `output` in place of `expected` to draw them from the results of a full run. The committed PDFs are included for comparison.

## Number check

The LaTeX manuscript sources are not included. If you have `main.tex` or `supplement.tex`, run:

```
python number_check/check_paper_numbers.py main.tex expected data
```

The script checks whether each printed number can match a value in a CSV or JSON file at the printed precision. It skips years and listed conventional values. A match does not establish that the number came from the correct row or analysis; use `manifest/analysis_manifest.csv` and the claim-status table for those links. Without the manuscript source, this check cannot be repeated from this repository alone.

The code is covered by the repository's MIT license.
