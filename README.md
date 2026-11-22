# NHIS 2023: Joint pain severity and recurrent falls

Analysis code for the manuscript *Joint pain severity and recurrent falls among US older adults reporting clinician-diagnosed arthritis or related conditions: a cross-sectional study of the 2023 National Health Interview Survey*.

## Data

Download the 2023 NHIS Sample Adult public-use CSV from the [official CDC/NCHS page](https://www.cdc.gov/nchs/nhis/documentation/2023-nhis.html) and extract `adult23.csv` into this repository directory. The dataset is deliberately excluded from this repository. No restricted-use NHIS data are required.

## Requirements

- Python 3.13.5 and the pinned dependencies in `requirements.txt` (these reproduce the manuscript's stated software versions; compatibility is not independently verified here).
- R and CRAN package `survey`. The R script exports its installed package information in `NHIS_R_survey_validation/R_sessionInfo.txt`.

## Commands

From the directory containing the scripts and CSV:

```bash
python -m pip install -r requirements.txt
python nhis2023_jointpain_falls_analysis.py
Rscript NHIS2023_R_survey_full_validation.R
python verify_outputs.py
```

The Python script writes results to `nhis_jointpain_project/`: sample-flow counts, population prevalence, Table 1 characteristics, Models 1–4, sensitivity/secondary outcomes, subgroup coefficients and interaction tests, VIF, three figures, a joint quadratic/cubic nonlinearity test and weighted Model 4-standardized predicted prevalences. The R script independently estimates survey-weighted prevalence and Models 1–4; provides fixed complete-case (N=4,296) results, inferential reference-distribution sensitivities and R/Python comparisons in `NHIS_R_survey_validation/`.

The script `verify_outputs.py` checks key selected results against the published manuscript values (including the newly implemented polynomial and standardized-prevalence computations). It will report `INCOMPLETE` if required files were not generated and `FAIL` if any selected estimates disagree with manuscript numbers. A PASS does **not** certify every supplementary-table cell or software environment.

## Statistical interpretation

The analysis uses the NHIS 2023 public-use sampling weight `WTFA_A`, strata `PSTRAT`, and clusters `PPSU` with domain-aware Taylor-linearized SEs. The nominal number 610 (662 PSUs minus 52 strata) is **not** an NCHS-certified effective inferential degrees-of-freedom value. The R script compares t(610), t(120), the normal reference and R's default residual degrees of freedom. Its embedded Python Model 4 reference numbers are historical comparison values, not inputs to the R model fits.

## Final R verification update (2026-11-23)

The R script was subsequently executed by the research team on Windows 11 with R 4.5.2 and `survey` 4.5. The five actual R output files are included in `NHIS_R_survey_validation/`. The R-versus-Python comparison file confirms Model 4 odds ratios agree to within 2.36e-10, with very small standard-error differences. Both fixed-N=4,296 regression and prevalence outputs have been supplied. The analyses were not rerun inside the packaging environment, so provenance of R execution is the user-supplied R output files and session information, rather than an independently executed R process here.

To check the provided outputs without rerunning R, run `python verify_outputs.py` after running Python with `adult23.csv` in the repository root. The outputs under `reference_python_outputs/` are provided reference outputs, not substitute raw source data.

## Prior Python audit status (2026-11-23)

**Python: EXECUTED on the user-supplied original `adult23.csv`.** The script completed with exit code 0 using the pinned versions listed above. The newly generated `python_reproducibility_audit.csv` reports **36/36 checks PASS**, including sample sizes, outcome counts, weighted prevalences, Model 4 ORs, sex/age interaction p values, polynomial Wald test, and standardized predicted prevalences. The source CSV SHA-256 was `441f2bc17761cb1c755b24c1af551e34b91e084da4cce8d78888a688b1acc5ab`. The data itself is **not** redistributed.

**R: NOT EXECUTED in the original Python packaging environment**, because `Rscript` was unavailable. Although the R script is included and explicitly implements independent model fits, the R/Python comparison remains to be run in an R environment. The user has now provided the R outputs; these are included in `NHIS_R_survey_validation/`.

Python-only recheck:

```bash
python nhis2023_jointpain_falls_analysis.py
python verify_python_outputs.py
```

Full independent R check:

```bash
Rscript NHIS2023_R_survey_full_validation.R
python verify_outputs.py
```

The audit does not certify every cell in supplementary tables S1–S7, independent inferential validity of the nominal 610 df convention, or portability to other software versions. Check the R session information and compare supplementary tables against calculated values before treating this as the final verified public release.
