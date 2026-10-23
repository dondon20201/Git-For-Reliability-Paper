# Analysis materials

Code for analysis and expected output from the paper, "How Many Problems Are Enough, and Which? A Reliability Ceiling for Individual-Difference Measures of Risky Choice" can be found in this repository. The folders contain the checks and associated materials:

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

The provenance file, `provenance/provenance.md`, lists dates from prior development. The specifications were added to the repository at a later time. The Git history of this repository therefore cannot verify the order of specifications and results listed in the provenance file. The specifications also contain references to the prior directory structure. For example, the source and render checks on the manuscript referenced in the integrity checks are not present in this repository.

## Figures

From the repository root, with Python 3, pandas, numpy, and matplotlib installed:

```
python figures/make_figures.py expected figures
```

The script uses the saved expected output and produces the figures `fig1_reliability.pdf`, `fig2_projection.pdf`, and `fig3_prior_studies.pdf`. To use the output from a complete run of the code rather than the saved expected output, replace `expected` with `output`. The saved PDFs are provided for reference.

## Number check

The LaTeX manuscript sources are not included. If you have `main.tex` or `supplement.tex`, run:

```
python number_check/check_paper_numbers.py main.tex expected data
```

The script will attempt to find a value in a JSON or CSV file that matches the printed value to the precision printed. It does not compare minus signs. Years and conventional values listed are ignored. Finding a match in a file does not imply that the value was produced in a specific analysis or row. Please see the analysis manifest in `manifest/analysis_manifest.csv` and the table of claims and their statuses for such information. This check cannot be reproduced using only this repository because the manuscript source is not available.

The code is covered by the repository's MIT license.
