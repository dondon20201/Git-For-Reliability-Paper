# Data

The raw data are not part of this repository. Put them in `data/raw/` with the
file names below. Steps 01, 09 and 12 download their datasets automatically if
the files are missing and stop if a checksum does not match. Most later steps
check the md5 of the CPC18 file again before reading it. The Frey et al. data
must be downloaded by hand; step 27 records the md5 of every Frey file it reads
in `output/elevation/results/p3_frey/provenance_manifest.csv`, so compare those
with the table below.

## 1. CPC18 estimation data (main dataset, steps 01-39)

Plonsky, O., Erev, I., & Ert, E. (2017). Calibration data for Choice
Prediction Competition 2018 (CPC18) (Version 1) [Data set]. Zenodo,
https://doi.org/10.5281/zenodo.845873 (record https://zenodo.org/records/845873),
CC-BY. Downloaded by step 01 if missing.

| File in `data/raw/` | md5 |
|---|---|
| `All estimation raw data.csv` | `38a26407a36c690be969f9d9d3fe3e5a` |
| `Dictionary for the choice prediction competition 2018 raw data.pdf` | `3980a6da798b4883167ebfcc7a2c52c5` |

`Overview of samples.pdf` from the same record (md5
`2a878c67521a87759c0567d1ea59b2b7`) is documentation only and is not read.

## 2. Official CPC18 baseline code (step 09)

Used to validate the reconstruction of multi-outcome lotteries. The official
"CPC18 BEASTsd R" files, as committed to https://github.com/naecker-lab/cpc18
at commit `b9caaca5874e3e388aea6580491ac639370b5bbd` (folder
`src/CPC18 BEASTsd R`), plus one cross-check file from
https://github.com/appgpk/choice_prediction at commit
`c2e7c571e7f0d7ca44c121265ae26a2d9df7e533` (`GetDistribution.py`, saved as
`GetDistribution_appgpk.py`). Downloaded by step 09 if missing. The exact URLs
are in `R/fetch_official_cpc18.R`.

| File in `data/raw/cpc18_official/` | md5 |
|---|---|
| `CPC18_getDist.R` | `338d7590949e11b598b5a8a081e4b79f` |
| `CPC18_get_pBetter.R` | `b7a70b8b15977cccfc71dc0aeb60eaba` |
| `distSample.R` | `389378a1febe93a1dbadced1a9fd57fa` |
| `CPC18_BEASTsd_simulation.R` | `c042c61d6149c5d939d4c9fd72a3bff9` |
| `CPC18_BEASTsd_pred.R` | `cd8d3430e9d473ee956c01173e550543` |
| `CPC18_EstSet210.csv` | `6c0817761e072ea4d642d7ff5909f839` |
| `CPC18_EstSet_fullDist.csv` | `f09bb9b51fcb2bb641a1b10a9ad6bad8` |
| `GetDistribution_appgpk.py` | `deab1ce2b5d2db8c18311201f396ecc3` |

## 3. choices13k (step 12, aggregate checks only)

Peterson, J. C., Bourgin, D. D., Agrawal, M., Reichman, D., & Griffiths, T. L.
(2021). Using large-scale experiments and machine learning to discover theories
of human decision-making. Science, 372, 1209-1214. Files from
https://github.com/jcpeterson/choices13k at commit
`821ae7e88386b508ebb46fae76fac63cb62ec876`. Downloaded by step 12 if missing.
These are problem-level aggregates and are used for aggregate agreement checks
only.

| File in `data/raw/choices13k/` | md5 |
|---|---|
| `c13k_selections.csv` | `65fd730812014d115c2cd0579f3122f9` |
| `c13k_problems.json` | `cc86e4e845dc2f48ff07c13b4663c20b` |
| `README.md` | `da9da2743206c7a1bc8517a84bb5614d` |

## 4. Frey et al. (2017) Basel-Berlin Risk Study (step 27, download by hand)

Frey, R., Pedroni, A., Mata, R., Rieskamp, J., & Hertwig, R. (2017). Risk
preference shares the psychometric structure of major psychological traits.
Science Advances, 3, e1701381. Data: OSF project https://osf.io/rce7g/.

Download the data and place the sample folders (`main/`, `retest_berlin/`,
`retest_basel/`) and `Codebook.pdf` in `data/raw/Frey Data/`. Step 27 reads the
files below; the other files of the download may stay in place. Step 27 also
checks with git that `data/raw/Frey Data` is ignored, so run it inside a git
clone of this repository.

| File in `data/raw/Frey Data/` | md5 |
|---|---|
| `main/measures/measures.Rdata` | `695d503790873261e38fd4697ccdad01` |
| `main/dfd/dfd_perpers.csv` | `37c7b7e584ae280a400c7b87635981a2` |
| `main/dfd/dfd_perprob.csv` | `b075ee28ec7e109fe25a71f75d9db4d0` |
| `main/dfe/dfe_perpers.csv` | `1b4f27a1dd5e2b18c260244bd6351f0c` |
| `main/dfe/dfe_perprob.csv` | `51e44af37330f144491c5686237bf34b` |
| `main/modelpars/modelpars.csv` | `a84911a8b71a5bfef6063af6b0ad2850` |
| `main/participants/participants.csv` | `e33da70b94437991f6158c4c82dcf9be` |
| `retest_berlin/measures/measures.Rdata` | `7a9222c2f1d7e8541668ccb5b67d72ce` |
| `retest_berlin/dfd/dfd_perpers.csv` | `26cfef773b588241171ac3960cd0c9fc` |
| `retest_berlin/dfd/dfd_perprob.csv` | `ef058058af9f9ab9a4063ca8d5838f75` |
| `retest_berlin/dfe/dfe_perpers.csv` | `326b39fb3140f3e21624d90f91ff711d` |
| `retest_berlin/dfe/dfe_perprob.csv` | `27e02614632dd246b19116c9f685bcf3` |
| `retest_berlin/modelpars/modelpars.csv` | `0b7f3e7bdb6d865c163cfd19a5da28ff` |
| `retest_berlin/participants/participants.csv` | `1a37b31249d9b8acc41ebd0acd93a29d` |
| `retest_basel/participants/participants.csv` | `3b9d0c74b008ca07c833a7a49e0357ad` |

## Input tables kept in this folder

- `prior_studies.csv`: published reliability results of prior studies (values,
  sources, design budgets), read by steps 08 and 13.
- `claim_prose.csv`: one-sentence interpretations joined to the claim-status
  table in step 14 (text only; no numbers).
