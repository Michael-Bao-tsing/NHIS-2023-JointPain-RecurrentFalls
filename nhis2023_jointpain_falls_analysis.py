#!/usr/bin/env python3
"""Reproducible analysis for 2023 NHIS joint pain severity and recurrent falls.

Study: U.S. adults aged >=65 years reporting clinician-diagnosed arthritis or related conditions.
Exposure: joint pain severity (0-3, 4-6, 7-10), with respondents reporting no
past-30-day joint symptoms assigned a pain score of 0, consistent with the
NHIS arthritis pain indicator construction.
Primary outcome: recurrent falls (>=2 falls in past 12 months).

Survey design: WTFA_A weights, PSTRAT strata, PPSU primary sampling units.
Taylor-series variance estimation is implemented by aggregating linearized
score/influence contributions within PSUs and strata. Domain estimates retain
all public-use strata/PSUs with zero contributions outside the analytic domain.

To run: Place adult23.csv in the same directory as this script, then execute:
    python nhis2023_jointpain_falls_analysis.py
Output files will be created in the nhis_jointpain_project subdirectory.
"""
from pathlib import Path
import json
import math
import warnings
import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm
from scipy import stats
from scipy.special import expit
from statsmodels.stats.outliers_influence import variance_inflation_factor
import matplotlib.pyplot as plt

# Portable input/output paths: independent of the current working directory.
SCRIPT_DIR = Path(__file__).resolve().parent
DATA = SCRIPT_DIR / 'adult23.csv'
OUT = SCRIPT_DIR / 'nhis_jointpain_project'

if not DATA.is_file():
    raise FileNotFoundError(
        f'NHIS source data not found: {DATA}. '
        'Download adult23.csv from the official 2023 NHIS Sample Adult public-use files '
        'and place it alongside this Python script.'
    )
OUT.mkdir(parents=True, exist_ok=True)

# ---------- Utility functions ----------
def yes_no(s):
    return s.map({1: 'Yes', 2: 'No'})

def survey_ratio_ci(full, domain, ycol, alpha=0.05):
    """Taylor-linearized weighted mean/proportion and t CI for a domain."""
    idx = domain & full[ycol].notna() & full['WTFA_A'].notna()
    w = full.loc[idx, 'WTFA_A'].astype(float)
    y = full.loc[idx, ycol].astype(float)
    est = float(np.sum(w*y)/np.sum(w))
    infl = pd.DataFrame({
        'PSTRAT': full['PSTRAT'], 'PPSU': full['PPSU'], 'u': 0.0
    })
    infl.loc[idx, 'u'] = w.values*(y.values-est)/float(np.sum(w))
    psu = infl.groupby(['PSTRAT','PPSU'], observed=True)['u'].sum().reset_index()
    var = 0.0
    for _, g in psu.groupby('PSTRAT', observed=True):
        vals = g['u'].to_numpy(float)
        m = len(vals)
        if m > 1:
            var += m/(m-1) * np.sum((vals-vals.mean())**2)
    npsu = full[['PSTRAT','PPSU']].drop_duplicates().shape[0]
    nstr = full['PSTRAT'].nunique()
    df = npsu - nstr
    se = math.sqrt(max(var,0))
    crit = stats.t.ppf(1-alpha/2, df)
    return est, se, max(0,est-crit*se), min(1,est+crit*se), df

def weighted_pct(data, mask):
    denom = data['WTFA_A'].sum()
    return 100*data.loc[mask, 'WTFA_A'].sum()/denom if denom > 0 else np.nan

def fit_survey_logit(full, domain, formula_rhs, outcome='recurrent'):
    """Survey-weighted logistic model with stratified PSU Taylor sandwich SE.

    Formula RHS is patsy syntax. Complete-case filtering is performed by patsy.
    Returns model point estimates, design-based covariance, OR table, and metadata.
    """
    base = full.loc[domain].copy()
    # Patsy drops rows with missing values.
    ymat, X = patsy.dmatrices(f'{outcome} ~ {formula_rhs}', base, return_type='dataframe', NA_action='drop')
    y = ymat.iloc[:,0].astype(float)
    d = base.loc[X.index].copy()
    w_raw = d['WTFA_A'].astype(float)
    w = w_raw / w_raw.mean()  # scaling leaves beta/sandwich invariant and improves numerics
    glm = sm.GLM(y, X, family=sm.families.Binomial(), freq_weights=w)
    res = glm.fit(maxiter=200, disp=0)
    beta = np.asarray(res.params)
    mu = expit(np.asarray(X) @ beta)
    Xv = np.asarray(X, float)
    wv = np.asarray(w, float)
    yv = np.asarray(y, float)
    # Observed Fisher information / pseudo-likelihood bread
    H = Xv.T @ ((wv*mu*(1-mu))[:,None] * Xv)
    Hinv = np.linalg.pinv(H)
    score_dom = (wv*(yv-mu))[:,None] * Xv

    # Domain estimation: preserve complete public-use design with zero score outside domain.
    p = Xv.shape[1]
    full_scores = pd.DataFrame(np.zeros((len(full),p)), index=full.index, columns=X.columns)
    full_scores.loc[X.index, :] = score_dom
    design = full[['PSTRAT','PPSU']].copy()
    design = pd.concat([design, full_scores], axis=1)
    psu = design.groupby(['PSTRAT','PPSU'], observed=True)[list(X.columns)].sum().reset_index()
    meat = np.zeros((p,p))
    for _, g in psu.groupby('PSTRAT', observed=True):
        U = g[X.columns].to_numpy(float)
        m = U.shape[0]
        if m > 1:
            C = U - U.mean(axis=0, keepdims=True)
            meat += m/(m-1) * (C.T @ C)
    cov = Hinv @ meat @ Hinv
    se = np.sqrt(np.clip(np.diag(cov),0,None))
    npsu = full[['PSTRAT','PPSU']].drop_duplicates().shape[0]
    nstr = full['PSTRAT'].nunique()
    ddf = npsu - nstr
    crit = stats.t.ppf(0.975, ddf)
    tval = beta/se
    pval = 2*stats.t.sf(np.abs(tval), ddf)
    out = pd.DataFrame({
        'term': X.columns,
        'beta': beta,
        'se': se,
        'OR': np.exp(beta),
        'CI_low': np.exp(beta-crit*se),
        'CI_high': np.exp(beta+crit*se),
        'p': pval,
    })
    meta = {
        'n': int(len(d)),
        'events': int(y.sum()),
        'weighted_population': float(w_raw.sum()),
        'design_df': int(ddf),
        'formula_rhs': formula_rhs,
        'outcome': outcome,
    }
    return out, cov, X, d, meta

def wald_joint(model_table, cov, terms, ddf):
    names = list(model_table['term'])
    inds = [names.index(t) for t in terms if t in names]
    if not inds:
        return np.nan, np.nan, 0
    b = model_table.set_index('term').loc[[names[i] for i in inds],'beta'].to_numpy()
    V = cov[np.ix_(inds,inds)]
    q = len(inds)
    W = float(b.T @ np.linalg.pinv(V) @ b)
    F = W/q
    p = float(stats.f.sf(F,q,ddf))
    return F,p,q

def fmt_or(row):
    return f"{row.OR:.2f} ({row.CI_low:.2f}-{row.CI_high:.2f})"

# ---------- Load and recode ----------
df = pd.read_csv(DATA, low_memory=False)
# design overview
nstr = df['PSTRAT'].nunique()
npsu = df[['PSTRAT','PPSU']].drop_duplicates().shape[0]
design_df = npsu-nstr

# Valid older age: codes 65..85, where 85 is 85+; 97/99 are nonresponses.
df['valid_age'] = df['AGEP_A'].between(18,85)
df['older'] = df['AGEP_A'].between(65,85)
df['arthritis'] = df['ARTHEV_A'].eq(1)
df['joint_symptom_valid'] = df['JNTSYMP_A'].isin([1,2])
# Assign 0 to those explicitly reporting no joint symptoms in past 30 days.
df['joint_pain'] = np.nan
df.loc[df['JNTSYMP_A'].eq(2),'joint_pain'] = 0.0
df.loc[df['JNTSYMP_A'].eq(1) & df['JNTPN_A'].between(0,10),'joint_pain'] = df.loc[df['JNTSYMP_A'].eq(1) & df['JNTPN_A'].between(0,10),'JNTPN_A']
df['pain_valid'] = df['joint_pain'].between(0,10)
df['pain_cat'] = pd.cut(df['joint_pain'], bins=[-0.01,3,6,10], labels=['0-3','4-6','7-10'])
df['pain_ord'] = df['pain_cat'].map({'0-3':0,'4-6':1,'7-10':2}).astype(float)

# Falls: BFALL12M_A 1 yes, 2 no. BFALLTIMES_A: 1 once, 2 twice, 3 3-4, 4 5-7, 5 8+.
df['fall_valid'] = df['BFALL12M_A'].eq(2) | (df['BFALL12M_A'].eq(1) & df['BFALLTIMES_A'].isin([1,2,3,4,5]))
df['any_fall'] = np.nan
df.loc[df['BFALL12M_A'].eq(2),'any_fall'] = 0
df.loc[df['BFALL12M_A'].eq(1) & df['BFALLTIMES_A'].isin([1,2,3,4,5]),'any_fall'] = 1
df['recurrent'] = np.nan
df.loc[df['BFALL12M_A'].eq(2),'recurrent'] = 0
df.loc[df['BFALL12M_A'].eq(1) & df['BFALLTIMES_A'].eq(1),'recurrent'] = 0
df.loc[df['BFALL12M_A'].eq(1) & df['BFALLTIMES_A'].isin([2,3,4,5]),'recurrent'] = 1

# Covariates
df['agegrp'] = pd.cut(df['AGEP_A'], bins=[64,74,84,85], labels=['65-74','75-84','85+'], include_lowest=True)
df['sex'] = df['SEX_A'].map({1:'Male',2:'Female'})

def raceeth(x):
    return {1:'Hispanic',2:'Non-Hispanic White',3:'Non-Hispanic Black',4:'Non-Hispanic Asian',5:'Other',6:'Other',7:'Other'}.get(x,np.nan)
df['raceeth'] = df['HISPALLP_A'].map(raceeth)

def education(x):
    if x in [1,2,3]: return '<High school/GED'
    if x==4: return 'High school graduate'
    if x in [5,6,7]: return 'Some college/associate'
    if x in [8,9,10]: return 'Bachelor or higher'
    return np.nan
df['education'] = df['EDUCP_A'].map(education)

def marital(x):
    if x in [1,2,3,8]: return 'Married/partnered'
    if x in [4,5,6,7]: return 'Not married/partnered'
    return np.nan
df['marital'] = df['MARSTAT_A'].map(marital)
df['bmi'] = df['BMICAT_A'].map({1:'Underweight',2:'Healthy weight',3:'Overweight',4:'Obesity'})
df['smoking'] = df['SMKCIGST_A'].map({1:'Current',2:'Current',3:'Former',4:'Never'})
df['phq2'] = df['PHQ2SCREEN_A'].map({1:'Positive',2:'Negative'})
df['vision_imp'] = df['VISIONDF_A'].map({1:'No',2:'Yes',3:'Yes',4:'Yes'})
df['hearing_imp'] = df['HEARINGDF_A'].map({1:'No',2:'Yes',3:'Yes',4:'Yes'})
for src,dst in [('HYPEV_A','hypertension'),('DIBEV_A','diabetes'),('COPDEV_A','copd'),('CHDEV_A','chd'),('STREV_A','stroke')]:
    df[dst] = df[src].map({1:'Yes',2:'No'})
df['balance'] = df['BALDIZZ_A'].map({1:'Yes',2:'No'})
df['mobility'] = df['DIFF_A'].map({1:'No difficulty',2:'Any difficulty',3:'Any difficulty',4:'Any difficulty'})
df['arth_limit'] = df['ARTHLMT_A'].map({1:'Yes',2:'No'})

# Primary domain
base0 = df['older']
base1 = base0 & df['arthritis']
base2 = base1 & df['joint_symptom_valid']
base3 = base2 & df['pain_valid']
primary = base3 & df['fall_valid']

flow = [
    ('2023 NHIS Sample Adult respondents', int(len(df)), None),
    ('Aged >=65 years (valid age codes 65-85)', int(base0.sum()), int(len(df)-base0.sum())),
    ('Doctor-diagnosed arthritis/related condition', int(base1.sum()), int(base0.sum()-base1.sum())),
    ('Valid past-30-day joint-symptom status', int(base2.sum()), int(base1.sum()-base2.sum())),
    ('Valid joint-pain exposure', int(base3.sum()), int(base2.sum()-base3.sum())),
    ('Valid fall outcome (primary analytic sample)', int(primary.sum()), int(base3.sum()-primary.sum())),
]
flow_df = pd.DataFrame(flow, columns=['Stage','Included_n','Excluded_from_previous_n'])
flow_df.to_csv(OUT/'flow_counts.csv', index=False)

# Verify published all-age prevalence construction.
all_arth = df['arthritis'] & df['joint_symptom_valid'] & df['pain_valid']
msjp = (df['joint_pain']>=4).astype(float)
verify = np.sum(df.loc[all_arth,'WTFA_A']*msjp.loc[all_arth]) / np.sum(df.loc[all_arth,'WTFA_A'])

# Primary prevalence
prev_rows=[]
for cat in ['0-3','4-6','7-10']:
    dom = primary & df['pain_cat'].eq(cat)
    est,se,lo,hi,_=survey_ratio_ci(df,dom,'recurrent')
    prev_rows.append({'Pain severity':cat,'Unweighted n':int(dom.sum()),'Recurrent falls n':int(df.loc[dom,'recurrent'].sum()),'Weighted prevalence %':100*est,'95% CI low %':100*lo,'95% CI high %':100*hi})
prev = pd.DataFrame(prev_rows)
prev.to_csv(OUT/'recurrent_falls_prevalence_by_pain.csv', index=False)

overall_rec = survey_ratio_ci(df,primary,'recurrent')
overall_any = survey_ratio_ci(df,primary,'any_fall')

# ---------- Table 1 ----------
table1_specs = [
    ('Age, years','agegrp',['65-74','75-84','85+']),
    ('Sex','sex',['Male','Female']),
    ('Race and ethnicity','raceeth',['Non-Hispanic White','Non-Hispanic Black','Hispanic','Non-Hispanic Asian','Other']),
    ('Education','education',['<High school/GED','High school graduate','Some college/associate','Bachelor or higher']),
    ('Marital status','marital',['Married/partnered','Not married/partnered']),
    ('BMI category','bmi',['Underweight','Healthy weight','Overweight','Obesity']),
    ('Smoking status','smoking',['Never','Former','Current']),
    ('PHQ-2 depression screen','phq2',['Negative','Positive']),
    ('Vision difficulty','vision_imp',['No','Yes']),
    ('Hearing difficulty','hearing_imp',['No','Yes']),
    ('Hypertension','hypertension',['No','Yes']),
    ('Diabetes','diabetes',['No','Yes']),
    ('COPD','copd',['No','Yes']),
    ('Coronary heart disease','chd',['No','Yes']),
    ('Stroke','stroke',['No','Yes']),
    ('Balance/dizziness problem, past 12 months','balance',['No','Yes']),
    ('Walking/climbing-steps difficulty','mobility',['No difficulty','Any difficulty']),
    ('Arthritis-attributable activity limitation','arth_limit',['No','Yes']),
]
cols = [('Overall',primary),('Pain 0-3',primary & df['pain_cat'].eq('0-3')),('Pain 4-6',primary & df['pain_cat'].eq('4-6')),('Pain 7-10',primary & df['pain_cat'].eq('7-10'))]
rows=[]
for label,var,levels in table1_specs:
    first=True
    for lev in levels:
        r={'Characteristic':label if first else '', 'Category':lev}
        first=False
        for cname,dom in cols:
            valid = dom & df[var].notna()
            n = int((valid & df[var].eq(lev)).sum())
            denom = df.loc[valid,'WTFA_A'].sum()
            pct = 100*df.loc[valid & df[var].eq(lev),'WTFA_A'].sum()/denom if denom else np.nan
            r[cname] = f'{n} ({pct:.1f})'
        rows.append(r)
t1=pd.DataFrame(rows)
t1.to_csv(OUT/'table1_characteristics.csv', index=False)

# ---------- Regression models ----------
# categorical reference ordering via explicit Treatment() formula
pain_term='C(pain_cat, Treatment(reference="0-3"))'
m1_rhs=pain_term
m2_rhs=pain_term+' + C(agegrp, Treatment(reference="65-74")) + C(sex, Treatment(reference="Male")) + C(raceeth, Treatment(reference="Non-Hispanic White")) + C(education, Treatment(reference="Bachelor or higher")) + C(marital, Treatment(reference="Married/partnered"))'
m3_rhs=m2_rhs+' + C(bmi, Treatment(reference="Healthy weight")) + C(smoking, Treatment(reference="Never")) + C(phq2, Treatment(reference="Negative")) + C(vision_imp, Treatment(reference="No")) + C(hearing_imp, Treatment(reference="No")) + C(hypertension, Treatment(reference="No")) + C(diabetes, Treatment(reference="No")) + C(copd, Treatment(reference="No")) + C(chd, Treatment(reference="No")) + C(stroke, Treatment(reference="No"))'
m4_rhs=m3_rhs+' + C(balance, Treatment(reference="No")) + C(mobility, Treatment(reference="No difficulty"))'
models={}
for nm,rhs in [('Model 1',m1_rhs),('Model 2',m2_rhs),('Model 3',m3_rhs),('Model 4',m4_rhs)]:
    tab,cov,X,d,meta=fit_survey_logit(df,primary,rhs)
    models[nm]=(tab,cov,X,d,meta)
    tab.to_csv(OUT/f'{nm.lower().replace(" ","_")}_full_coefficients.csv', index=False)

# Main pain-only rows across models.
main_rows=[]
for nm,(tab,cov,X,d,meta) in models.items():
    for cat in ['4-6','7-10']:
        key=f'{pain_term}[T.{cat}]'
        rr=tab[tab.term.eq(key)].iloc[0]
        main_rows.append({'Model':nm,'Pain severity':cat,'OR':rr.OR,'CI_low':rr.CI_low,'CI_high':rr.CI_high,'p':rr.p,'N':meta['n'],'Events':meta['events']})
main=pd.DataFrame(main_rows)
main.to_csv(OUT/'table2_main_models.csv', index=False)

# Trend (ordered categories) and continuous score, fully adjusted.
trend_rhs=m4_rhs.replace(pain_term,'pain_ord')
trend_tab,trend_cov,trend_X,trend_d,trend_meta=fit_survey_logit(df,primary,trend_rhs)
cont_rhs=m4_rhs.replace(pain_term,'joint_pain')
cont_tab,cont_cov,cont_X,cont_d,cont_meta=fit_survey_logit(df,primary,cont_rhs)
trend_row=trend_tab[trend_tab.term.eq('pain_ord')].iloc[0]
cont_row=cont_tab[cont_tab.term.eq('joint_pain')].iloc[0]

# Sensitivity 1: symptomatic only
symptom_domain=primary & df['JNTSYMP_A'].eq(1)
sym_tab,sym_cov,sym_X,sym_d,sym_meta=fit_survey_logit(df,symptom_domain,m4_rhs)
# Sensitivity 2: further adjustment for arthritis-attributable limitation
lim_rhs=m4_rhs+' + C(arth_limit, Treatment(reference="No"))'
lim_tab,lim_cov,lim_X,lim_d,lim_meta=fit_survey_logit(df,primary,lim_rhs)
# Secondary: any fall
any_domain=base3 & df['any_fall'].notna()
any_tab,any_cov,any_X,any_d,any_meta=fit_survey_logit(df,any_domain,m4_rhs,outcome='any_fall')

# Collect Table 3 sensitivity/secondary
sens=[]
def add_sens(label, tab, meta, term, unit):
    rr=tab[tab.term.eq(term)].iloc[0]
    sens.append({'Analysis':label,'Exposure contrast':unit,'OR':rr.OR,'CI_low':rr.CI_low,'CI_high':rr.CI_high,'p':rr.p,'N':meta['n'],'Events':meta['events']})
add_sens('Linear trend across pain categories',trend_tab,trend_meta,'pain_ord','Per one-category increase')
add_sens('Continuous pain score',cont_tab,cont_meta,'joint_pain','Per 1-point increase')
add_sens('Exclude respondents with no joint symptoms',sym_tab,sym_meta,f'{pain_term}[T.7-10]','Severe vs no/mild')
add_sens('Additional adjustment for arthritis-attributable activity limitation',lim_tab,lim_meta,f'{pain_term}[T.7-10]','Severe vs no/mild')
add_sens('Secondary outcome: any fall',any_tab,any_meta,f'{pain_term}[T.7-10]','Severe vs no/mild')
sens_df=pd.DataFrame(sens)
sens_df.to_csv(OUT/'table3_sensitivity_analyses.csv',index=False)

# Full Model 4 supplementary table
models['Model 4'][0].to_csv(OUT/'supp_table_full_model4.csv', index=False)

# Subgroups: sex and age; remove stratifier from covariate set.
sex_rhs=pain_term+' + C(agegrp, Treatment(reference="65-74")) + C(raceeth, Treatment(reference="Non-Hispanic White")) + C(education, Treatment(reference="Bachelor or higher")) + C(marital, Treatment(reference="Married/partnered")) + C(bmi, Treatment(reference="Healthy weight")) + C(smoking, Treatment(reference="Never")) + C(phq2, Treatment(reference="Negative")) + C(vision_imp, Treatment(reference="No")) + C(hearing_imp, Treatment(reference="No")) + C(hypertension, Treatment(reference="No")) + C(diabetes, Treatment(reference="No")) + C(copd, Treatment(reference="No")) + C(chd, Treatment(reference="No")) + C(stroke, Treatment(reference="No")) + C(balance, Treatment(reference="No")) + C(mobility, Treatment(reference="No difficulty"))'
age_rhs=pain_term+' + C(sex, Treatment(reference="Male")) + C(raceeth, Treatment(reference="Non-Hispanic White")) + C(education, Treatment(reference="Bachelor or higher")) + C(marital, Treatment(reference="Married/partnered")) + C(bmi, Treatment(reference="Healthy weight")) + C(smoking, Treatment(reference="Never")) + C(phq2, Treatment(reference="Negative")) + C(vision_imp, Treatment(reference="No")) + C(hearing_imp, Treatment(reference="No")) + C(hypertension, Treatment(reference="No")) + C(diabetes, Treatment(reference="No")) + C(copd, Treatment(reference="No")) + C(chd, Treatment(reference="No")) + C(stroke, Treatment(reference="No")) + C(balance, Treatment(reference="No")) + C(mobility, Treatment(reference="No difficulty"))'
sub=[]
for val in ['Male','Female']:
    tab,cov,X,d,meta=fit_survey_logit(df,primary & df['sex'].eq(val),sex_rhs)
    rr=tab[tab.term.eq(f'{pain_term}[T.7-10]')].iloc[0]
    sub.append({'Stratifier':'Sex','Subgroup':val,'N':meta['n'],'Events':meta['events'],'OR':rr.OR,'CI_low':rr.CI_low,'CI_high':rr.CI_high,'p':rr.p})
# age 65-74 vs 75+
df['age_2group'] = pd.Series(pd.NA, index=df.index, dtype='object')
df.loc[df['AGEP_A'].between(65,74),'age_2group']='65-74'
df.loc[df['AGEP_A'].between(75,85),'age_2group']='75+'
for val in ['65-74','75+']:
    tab,cov,X,d,meta=fit_survey_logit(df,primary & df['age_2group'].eq(val),age_rhs)
    rr=tab[tab.term.eq(f'{pain_term}[T.7-10]')].iloc[0]
    sub.append({'Stratifier':'Age','Subgroup':val,'N':meta['n'],'Events':meta['events'],'OR':rr.OR,'CI_low':rr.CI_low,'CI_high':rr.CI_high,'p':rr.p})
sub_df=pd.DataFrame(sub)

# Interaction models and joint tests.
# Sex interaction (model replaces pain+sex with their interaction; all other full covariates).
sex_int_rhs = pain_term+' * C(sex, Treatment(reference=\"Male\")) + C(agegrp, Treatment(reference=\"65-74\")) + C(raceeth, Treatment(reference=\"Non-Hispanic White\")) + C(education, Treatment(reference=\"Bachelor or higher\")) + C(marital, Treatment(reference=\"Married/partnered\")) + C(bmi, Treatment(reference=\"Healthy weight\")) + C(smoking, Treatment(reference=\"Never\")) + C(phq2, Treatment(reference=\"Negative\")) + C(vision_imp, Treatment(reference=\"No\")) + C(hearing_imp, Treatment(reference=\"No\")) + C(hypertension, Treatment(reference=\"No\")) + C(diabetes, Treatment(reference=\"No\")) + C(copd, Treatment(reference=\"No\")) + C(chd, Treatment(reference=\"No\")) + C(stroke, Treatment(reference=\"No\")) + C(balance, Treatment(reference=\"No\")) + C(mobility, Treatment(reference=\"No difficulty\"))'
sex_int_tab,sex_int_cov,sex_int_X,sex_int_d,sex_int_meta=fit_survey_logit(df,primary,sex_int_rhs)
sex_inter_terms=[t for t in sex_int_tab.term if ':' in t and 'pain_cat' in t and 'sex' in t]
_,sex_int_p,_=wald_joint(sex_int_tab,sex_int_cov,sex_inter_terms,sex_int_meta['design_df'])
# age interaction using two-group age; create a full formula without agegrp and with age_2group.
base_no_age=m4_rhs.replace(' + C(agegrp, Treatment(reference="65-74"))','')
age_int_rhs=base_no_age.replace(pain_term, f'{pain_term} * C(age_2group, Treatment(reference="65-74"))')
age_int_tab,age_int_cov,age_int_X,age_int_d,age_int_meta=fit_survey_logit(df,primary,age_int_rhs)
age_inter_terms=[t for t in age_int_tab.term if ':' in t and 'pain_cat' in t and 'age_2group' in t]
_,age_int_p,_=wald_joint(age_int_tab,age_int_cov,age_inter_terms,age_int_meta['design_df'])
sub_df['Interaction p']=np.nan
sub_df.loc[(sub_df.Stratifier=='Sex') & (sub_df.groupby('Stratifier').cumcount()==0),'Interaction p']=sex_int_p
sub_df.loc[(sub_df.Stratifier=='Age') & (sub_df.groupby('Stratifier').cumcount()==0),'Interaction p']=age_int_p
sub_df.to_csv(OUT/'subgroup_analyses.csv',index=False)

# VIF for Model 4 matrix (exclude intercept).
X4=models['Model 4'][2]
# VIF may be high for dummy sets but useful QC.
vif=[]
for i,c in enumerate(X4.columns):
    if c=='Intercept': continue
    try: v=variance_inflation_factor(np.asarray(X4,float),i)
    except Exception: v=np.nan
    vif.append({'term':c,'VIF':v})
vif_df=pd.DataFrame(vif)
vif_df.to_csv(OUT/'vif_model4.csv',index=False)

# ---------- Figures ----------
plt.rcParams.update({'font.size':10,'font.family':'DejaVu Sans'})
# Figure 1: participant flow diagram
fig,ax=plt.subplots(figsize=(7.2,8.8))
ax.axis('off')
ys=np.linspace(0.90,0.10,len(flow_df))
for i,(y,(_,r)) in enumerate(zip(ys,flow_df.iterrows())):
    txt=f"{r['Stage']}\nN = {r['Included_n']:,}"
    if i>0 and pd.notna(r['Excluded_from_previous_n']):
        txt += f"\nExcluded from previous step: {int(r['Excluded_from_previous_n']):,}"
    ax.text(0.5,y,txt,ha='center',va='center',bbox=dict(boxstyle='round,pad=0.55',fc='white',ec='black',lw=1.0),transform=ax.transAxes)
    if i<len(ys)-1:
        ax.annotate('',xy=(0.5,ys[i+1]+0.055),xytext=(0.5,y-0.055),xycoords=ax.transAxes,textcoords=ax.transAxes,arrowprops=dict(arrowstyle='->',lw=1.0))
ax.set_title('Figure 1. Participant selection from the 2023 National Health Interview Survey',fontsize=12,pad=12)
fig.text(0.5,0.015,'Data source: National Center for Health Statistics, National Health Interview Survey, 2023.',ha='center',fontsize=8)
fig.tight_layout(rect=[0,0.035,1,1])
fig.savefig(OUT/'Figure1_participant_flow.png',dpi=300,bbox_inches='tight')
plt.close(fig)

# Figure 2: recurrent fall prevalence by pain severity
fig,ax=plt.subplots(figsize=(6.6,4.6))
x=np.arange(3)
y=prev['Weighted prevalence %'].to_numpy()
lo=prev['95% CI low %'].to_numpy(); hi=prev['95% CI high %'].to_numpy()
ax.errorbar(x,y,yerr=[y-lo,hi-y],fmt='o',capsize=5,markersize=7,linewidth=1.5)
ax.set_xticks(x, ['0-3\nNo/mild','4-6\nModerate','7-10\nSevere'])
ax.set_ylabel('Survey-weighted prevalence of recurrent falls (%)')
ax.set_xlabel('Average joint pain severity, past 30 days')
ax.set_ylim(0,max(30,hi.max()+3))
ax.grid(axis='y',alpha=0.25)
for xi,yi in zip(x,y): ax.text(xi,yi+1.6,f'{yi:.1f}%',ha='center',va='bottom')
fig.text(0.5,0.01,'Data source: National Center for Health Statistics, National Health Interview Survey, 2023.',ha='center',fontsize=8)
fig.tight_layout(rect=[0,0.04,1,1])
fig.savefig(OUT/'Figure2_recurrent_falls_prevalence.png',dpi=300,bbox_inches='tight')
plt.close(fig)

# Figure 3: fully adjusted and sensitivity results severe vs no/mild
forest=[
    ('Primary: recurrent falls, Model 4', main[(main.Model=='Model 4') & (main['Pain severity']=='7-10')].iloc[0]),
    ('Exclude no joint symptoms', sens_df[sens_df.Analysis.str.startswith('Exclude')].iloc[0]),
    ('+ Arthritis activity limitation', sens_df[sens_df.Analysis.str.startswith('Additional')].iloc[0]),
    ('Secondary outcome: any fall', sens_df[sens_df.Analysis.str.startswith('Secondary')].iloc[0]),
]
labels=[a for a,b in forest]
ors=np.array([float(b.OR) for a,b in forest]); los=np.array([float(b.CI_low) for a,b in forest]); his=np.array([float(b.CI_high) for a,b in forest])
ypos=np.arange(len(labels))[::-1]
fig,ax=plt.subplots(figsize=(7.2,4.3))
ax.errorbar(ors,ypos,xerr=[ors-los,his-ors],fmt='o',capsize=4,linewidth=1.4)
ax.axvline(1,color='black',lw=1,ls='--')
ax.set_yticks(ypos,labels)
ax.set_xlabel('Odds ratio for severe joint pain (7-10 vs 0-3)')
ax.set_xlim(0.75,max(2.6,his.max()+0.15))
ax.grid(axis='x',alpha=0.2)
for o,l,h,yv in zip(ors,los,his,ypos):
    ax.text(h+0.04,yv,f'{o:.2f} ({l:.2f}-{h:.2f})',va='center',fontsize=9)
fig.text(0.5,0.01,'Data source: National Center for Health Statistics, National Health Interview Survey, 2023.',ha='center',fontsize=8)
fig.tight_layout(rect=[0,0.04,1,1])
fig.savefig(OUT/'Figure3_sensitivity_forest.png',dpi=300,bbox_inches='tight')
plt.close(fig)

# ---------- Manuscript additional analyses (Supplementary Table S7) ----------
# Same Model 4 complete-case sample, with continuous pain score instead of categories.
# The joint test of quadratic and cubic terms probes departure from linear log-odds.
df['pain_squared'] = df['joint_pain'] ** 2
df['pain_cubed'] = df['joint_pain'] ** 3
poly_rhs = cont_rhs + ' + pain_squared + pain_cubed'
poly_tab, poly_cov, poly_X, poly_d, poly_meta = fit_survey_logit(df, primary, poly_rhs)
poly_F, poly_p, poly_q = wald_joint(poly_tab, poly_cov,
                                   ['pain_squared', 'pain_cubed'],
                                   poly_meta['design_df'])
if poly_meta['n'] != cont_meta['n']:
    raise AssertionError('Polynomial and linear pain models must use the same sample')
nonlinear_summary = pd.DataFrame([{
    'test': 'Quadratic and cubic pain terms (joint design-based Wald F)',
    'F': poly_F, 'numerator_df': poly_q,
    'denominator_df_nominal': poly_meta['design_df'],
    'p': poly_p, 'N': poly_meta['n']
}])
nonlinear_summary.to_csv(OUT/'supp_nonlinearity_test.csv', index=False)

# Weighted model-standardized predictions in the Model 4 complete-case cohort.
# Use the same formula design_info and model coefficients as the fitted Model 4;
# predict for each respondent while setting their exposure category to the target.
model4_tab, model4_cov, model4_X, model4_d, model4_meta = models['Model 4']
standardized_rows = []
for category in ['0-3', '4-6', '7-10']:
    counterfactual = model4_d.copy()
    counterfactual['pain_cat'] = pd.Categorical(
        [category] * len(counterfactual),
        categories=['0-3', '4-6', '7-10'], ordered=True)
    X_new = patsy.build_design_matrices(
        [model4_X.design_info], counterfactual, return_type='dataframe')[0]
    X_new = X_new[model4_X.columns]
    predicted = expit(X_new.to_numpy(float) @ model4_tab['beta'].to_numpy(float))
    average = np.average(predicted, weights=model4_d['WTFA_A'].to_numpy(float))
    standardized_rows.append({
        'Pain severity': category, 'N': model4_meta['n'],
        'Model-standardized predicted recurrent fall prevalence %': 100*float(average)
    })
standardized_df = pd.DataFrame(standardized_rows)
standardized_df.to_csv(OUT/'supp_model_standardized_prevalence.csv', index=False)

# ---------- Summary JSON / text ----------
summary={
    'source_n':int(len(df)),
    'older_n':int(base0.sum()),
    'arthritis_older_n':int(base1.sum()),
    'primary_n':int(primary.sum()),
    'primary_recurrent_events':int(df.loc[primary,'recurrent'].sum()),
    'primary_any_fall_events':int(df.loc[primary,'any_fall'].sum()),
    'weighted_population':float(df.loc[primary,'WTFA_A'].sum()),
    'weighted_recurrent_prevalence':float(overall_rec[0]),
    'weighted_recurrent_ci':[float(overall_rec[2]),float(overall_rec[3])],
    'weighted_any_fall_prevalence':float(overall_any[0]),
    'weighted_any_fall_ci':[float(overall_any[2]),float(overall_any[3])],
    'design_strata':int(nstr),'design_psu':int(npsu),'design_df':int(design_df),
    'replication_check_all_age_n':int(all_arth.sum()),
    'replication_check_ms_joint_pain_pct':float(100*verify),
    'model4_n':models['Model 4'][4]['n'],
    'model4_events':models['Model 4'][4]['events'],
    'model4_moderate':main[(main.Model=='Model 4')&(main['Pain severity']=='4-6')].iloc[0][['OR','CI_low','CI_high','p']].to_dict(),
    'model4_severe':main[(main.Model=='Model 4')&(main['Pain severity']=='7-10')].iloc[0][['OR','CI_low','CI_high','p']].to_dict(),
    'trend':trend_row[['OR','CI_low','CI_high','p']].to_dict(),
    'continuous':cont_row[['OR','CI_low','CI_high','p']].to_dict(),
    'sex_interaction_p':sex_int_p,'age_interaction_p':age_int_p,
    'max_vif':float(vif_df.VIF.max()),
    'nonlinearity_F':float(poly_F),
    'nonlinearity_p':float(poly_p),
    'model_standardized_prevalence_pct':standardized_df['Model-standardized predicted recurrent fall prevalence %'].to_list(),
}
with open(OUT/'analysis_summary.json','w') as f: json.dump(summary,f,indent=2)

print(json.dumps(summary,indent=2))
print('\nPrevalence by pain:\n',prev.to_string(index=False))
print('\nMain models:\n',main.to_string(index=False))
print('\nSensitivity:\n',sens_df.to_string(index=False))
print('\nSubgroups:\n',sub_df.to_string(index=False))
