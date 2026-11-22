#!/usr/bin/env python3
"""Numerical audit of independently regenerated Python results against manuscript."""
from pathlib import Path
import json, math, csv, sys
import pandas as pd
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'nhis_jointpain_project'
s=json.loads((OUT/'analysis_summary.json').read_text())
rows=[]
def test(label,actual,expected,tolerance):
    ok=math.isfinite(float(actual)) and abs(float(actual)-expected)<=tolerance
    rows.append((label,'PASS' if ok else 'FAIL',actual,expected,tolerance))
for name,key,expect in [('Full NHIS sample','source_n',29522),('Older with arthritis','arthritis_older_n',4654),('Eligible sample','primary_n',4473),('Recurrent events','primary_recurrent_events',709),('Any fall events','primary_any_fall_events',1387),('Model 4 N','model4_n',4296),('Strata','design_strata',52),('PSUs','design_psu',662),('Nominal design df','design_df',610)]:test(name,s[key],expect,0)
test('Weighted target total',s['weighted_population'],26510428,0.5)
test('Overall recurrent prevalence %',100*s['weighted_recurrent_prevalence'],15.0,.05)
test('Overall any fall prevalence %',100*s['weighted_any_fall_prevalence'],30.5,.05)
test('Model 4 moderate OR',s['model4_moderate']['OR'],1.2735382805821989,1e-6)
test('Model 4 severe OR',s['model4_severe']['OR'],1.747145266960071,1e-6)
test('Ordinal OR',s['trend']['OR'],1.32,.005)
test('Continuous OR',s['continuous']['OR'],1.08,.005)
test('Sex interaction p',s['sex_interaction_p'],.286,.0005)
test('Age interaction p',s['age_interaction_p'],.986,.0005)
test('VIF max',s['max_vif'],1.74,.005)
test('Nonlinearity F',s['nonlinearity_F'],.140,.0005)
test('Nonlinearity p',s['nonlinearity_p'],.870,.0005)
p=pd.read_csv(OUT/'recurrent_falls_prevalence_by_pain.csv')
for cat,expect,n,events in [('0-3',10.2,1996,216),('4-6',15.3,1464,240),('7-10',23.5,1013,253)]:
    r=p.loc[p['Pain severity'].astype(str)==cat].iloc[0]
    test(cat+' weighted prevalence %',r['Weighted prevalence %'],expect,.05)
    test(cat+' N',r['Unweighted n'],n,0)
    test(cat+' events',r['Recurrent falls n'],events,0)
for cat,actual,expect in zip(['0-3','4-6','7-10'],s['model_standardized_prevalence_pct'],[12.4,14.9,18.7]):test(cat+' standardized prevalence',actual,expect,.05)
a=pd.read_csv(OUT/'table3_sensitivity_analyses.csv')
for prefix,expect in [('Exclude respondents',1.67),('Additional adjustment',1.64),('Secondary outcome',1.31)]:
    r=a[a['Analysis'].str.startswith(prefix)].iloc[0]
    test(prefix+' severe OR',r.OR,expect,.005)
report=pd.DataFrame(rows,columns=['check','status','observed','manuscript','tolerance'])
report.to_csv(ROOT/'python_reproducibility_audit.csv',index=False)
print(report.to_string(index=False))
print(f'\n{sum(report.status=="PASS")}/{len(report)} checks passed')
sys.exit(0 if (report.status=='PASS').all() else 1)
