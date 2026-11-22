#!/usr/bin/env python3
"""Fail-loud audit of key manuscript numbers against reproducibly generated CSV/JSON.
Run after Python and R analyses. This does not replace line-by-line review of tables.
"""
from pathlib import Path
import json
import math
import sys
import pandas as pd
ROOT = Path(__file__).resolve().parent
PYOUT = ROOT/'nhis_jointpain_project'
ROUT = ROOT/'NHIS_R_survey_validation'
checks=[]
def check(name, value, expected, tol):
    ok = math.isfinite(float(value)) and abs(float(value)-expected) <= tol
    checks.append((name, 'PASS' if ok else 'FAIL',float(value),expected,tol))

def load_csv(dir, name):
    f=dir/name
    if not f.is_file():
        raise FileNotFoundError(f'Missing {f}; run analysis scripts first')
    return pd.read_csv(f)

try:
    summary = json.loads((PYOUT/'analysis_summary.json').read_text())
    for name,k,val in [('full sample','source_n',29522),('eligible sample','primary_n',4473),
                       ('recurrent events','primary_recurrent_events',709),('Model 4 complete cases','model4_n',4296),
                       ('design strata','design_strata',52),('stratum–PSUs','design_psu',662)]:
        check(name,summary[k],val,0)
    check('severe pain Model 4 OR',summary['model4_severe']['OR'],1.747145266960071,0.000005)
    check('moderate pain Model 4 OR',summary['model4_moderate']['OR'],1.2735382805821989,0.000005)
    check('overall recurrent prevalence %',100*summary['weighted_recurrent_prevalence'],15.0,0.05)
    prev=load_csv(PYOUT,'recurrent_falls_prevalence_by_pain.csv')
    for cat,target in [('0-3',10.2),('4-6',15.3),('7-10',23.5)]:
        row=prev.loc[prev['Pain severity'].astype(str).eq(cat)].iloc[0]
        check(f'prevalence % {cat}',row['Weighted prevalence %'],target,0.05)
    nonlinear=load_csv(PYOUT,'supp_nonlinearity_test.csv').iloc[0]
    check('nonlinearity F (rounded)',nonlinear.F,0.140,0.0005)
    check('nonlinearity p (rounded)',nonlinear.p,0.870,0.0005)
    std=load_csv(PYOUT,'supp_model_standardized_prevalence.csv')
    for cat,target in [('0-3',12.4),('4-6',14.9),('7-10',18.7)]:
        row=std.loc[std['Pain severity'].astype(str).eq(cat)].iloc[0]
        check(f'standardized prevalence % {cat}',row['Model-standardized predicted recurrent fall prevalence %'],target,0.05)
    r=load_csv(ROUT,'R_vs_Python_Model4_comparison.csv')
    for _,row in r.iterrows():
        check(f'R/Python Model 4 OR agreement {row.contrast}',row['absolute_OR_difference'],0,1e-5)
        check(f'R/Python Model 4 SE relative difference (%) {row.contrast}',row['relative_SE_difference_pct'],0,0.01)
    fixed=load_csv(ROUT,'R_fixed_4296_models_all_df_options.csv')
    check('fixed sample rows',len(fixed),32,0)  # 4 models × 2 pain contrasts × 4 inference rules
    report=pd.DataFrame(checks,columns=['check','status','observed','expected','tolerance'])
    report.to_csv(ROOT/'reproducibility_check_results.csv',index=False)
    print(report.to_string(index=False))
    if any(x[1]=='FAIL' for x in checks): sys.exit('At least one manuscript value failed reproducibility verification')
    print('\nAll selected numerical checks PASS. Full table-level audit is still required.')
except (FileNotFoundError,KeyError,IndexError) as e:
    sys.exit(f'INCOMPLETE: {e}')
