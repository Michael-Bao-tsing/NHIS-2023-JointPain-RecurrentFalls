# NHIS 2023 reproducibility review — 2026-11-23

## Executed tests
- Source data: original user-supplied `adult23.csv`, 29,522 rows.
- SHA-256: `441f2bc17761cb1c755b24c1af551e34b91e084da4cce8d78888a688b1acc5ab`.
- Python executable: completed successfully (exit code 0), all 36 prespecified selected-result checks passed (see `python_reproducibility_audit.csv`).
- Python dependencies used: pandas 2.2.3; numpy 2.3.5; scipy 1.17.0; statsmodels 0.14.6; patsy 1.0.2; matplotlib 3.10.8.

## Confirmed results
- Primary analytic N = 4473, recurrent-fall events = 709, full Model 4 N = 4296.
- Model 4 OR moderate = 1.2735382806; severe = 1.7471452670.
- Nonlinearity joint Wald F = 0.13956171035; p = 0.86976711262 (rounded 0.140; 0.870).
- Model-standardized weighted recurrent-fall prevalence = 12.3686176%, 14.8827883%, 18.7223469% (rounded 12.4%, 14.9%, 18.7%).

## Outstanding
- R script cannot be executed without R + `survey`; no new R/Python concordance claim is made here.
- The combined `verify_outputs.py` returns INCOMPLETE until the R outputs exist.
- Complete line-by-line comparison with formatted supplementary tables has not been performed.
- Code does not establish an NHIS-certified effective degrees of freedom; the t(610) convention remains a sensitivity/interpretation caveat.

## Sharing
- The original NHIS microdata CSV is excluded from this code package; obtain from CDC/NCHS.
- No public GitHub URL has been created or verified.

## Follow-up R execution
The research team provided five R outputs on 2026-11-23 from R 4.5.2 / survey 4.5 on Windows 11. Outputs are included in NHIS_R_survey_validation/. Model 4 R-versus-Python maximum absolute OR difference is 2.350893e-10. These files have been inspected; the packaging environment did not itself execute R.
