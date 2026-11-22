#!/usr/bin/env Rscript
# ============================================================================
# NHIS 2023: independent R/survey validation of the joint-pain/falls paper
# Corresponds to: nhis2023_jointpain_falls_analysis.py and adult23.csv
# Outputs: independent weighted prevalence, Model 1-4 OR/SE/CI/p,
#          identical-sample (N=4,296) models, Python comparisons, df sensitivity.
# R package required: survey. No R model results are hard-coded.
# IMPORTANT: 610 = public-use PSU-strata bookkeeping df, not an NCHS-certified
#            effective df for all NHIS outcomes. 120 is sensitivity ONLY.
# ============================================================================

# 1. SETUP --------------------------------------------------------------------
if (!requireNamespace("survey", quietly = TRUE)) {
  stop("Package 'survey' missing. Install first: install.packages('survey')")
}
suppressPackageStartupMessages(library(survey))
options(survey.lonely.psu = "adjust")
# For domain estimation, subsetting survey design (not raw data) is essential.
# This option controls handling if a domain contains a singleton PSU.
options(survey.adjust.domain.lonely = TRUE)

args_all <- commandArgs(trailingOnly=FALSE)
file_arg <- grep('^--file=', args_all, value=TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub('^--file=', '', file_arg[1]))) else getwd()
files <- file.path(script_dir, c("adult23.csv", "adult23(1).csv"))
infile <- files[file.exists(files)][1]
if (is.na(infile)) stop("Place adult23.csv alongside this R script, then retry.")
out_dir <- file.path(script_dir, "NHIS_R_survey_validation")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat("Input:", normalizePath(infile), "\n")
d <- read.csv(infile, check.names = FALSE, stringsAsFactors = FALSE)
needed <- c("WTFA_A","PSTRAT","PPSU","AGEP_A","ARTHEV_A",
            "JNTSYMP_A","JNTPN_A","BFALL12M_A","BFALLTIMES_A",
            "SEX_A","HISPALLP_A","EDUCP_A","MARSTAT_A","BMICAT_A",
            "SMKCIGST_A","PHQ2SCREEN_A","VISIONDF_A","HEARINGDF_A",
            "HYPEV_A","DIBEV_A","COPDEV_A","CHDEV_A","STREV_A",
            "BALDIZZ_A","DIFF_A")
if (!all(needed %in% names(d))) {
  stop("Missing NHIS columns: ", paste(setdiff(needed, names(d)), collapse=", "))
}
stopifnot(!anyNA(d[c("WTFA_A","PSTRAT","PPSU")]), all(d$WTFA_A > 0))

# 2. RECODE: EXACT SAME RULES AS THE SUPPLIED PYTHON SCRIPT --------------------
# A factor mapping helper; unspecified codes become NA (not a valid category).
fmap <- function(x, codes, levels) {
  factor(unname(codes[as.character(x)]), levels = levels)
}

d$older <- d$AGEP_A %in% 65:85
d$arthritis <- d$ARTHEV_A %in% 1
valid_symptom <- d$JNTSYMP_A %in% c(1,2)

d$joint_pain <- rep(NA_real_, nrow(d))
d$joint_pain[d$JNTSYMP_A %in% 2] <- 0
ix <- d$JNTSYMP_A %in% 1 & d$JNTPN_A >= 0 & d$JNTPN_A <= 10
ix[is.na(ix)] <- FALSE
d$joint_pain[ix] <- d$JNTPN_A[ix]
d$pain_cat <- cut(d$joint_pain, breaks=c(-0.01,3,6,10),
                  labels=c("0-3","4-6","7-10"), include.lowest=TRUE)
d$pain_cat <- factor(d$pain_cat, levels=c("0-3","4-6","7-10"))

# BFALLTIMES_A is categorical: 1=one, 2=two, 3=3-4, 4=5-7, 5=8+
# Codes 7,8,9 are not valid fall counts and must not be >=2 events.
fall_valid <- d$BFALL12M_A %in% 2 |
  (d$BFALL12M_A %in% 1 & d$BFALLTIMES_A %in% 1:5)
d$recurrent <- rep(NA_real_, nrow(d))
d$recurrent[d$BFALL12M_A %in% 2] <- 0
d$recurrent[d$BFALL12M_A %in% 1 & d$BFALLTIMES_A %in% 1] <- 0
d$recurrent[d$BFALL12M_A %in% 1 & d$BFALLTIMES_A %in% 2:5] <- 1

d$primary <- as.integer(d$older & d$arthritis & valid_symptom &
                        !is.na(d$joint_pain) & fall_valid)

# Reference levels exactly reproduce the Python patsy Treatment() references.
d$agegrp <- factor(ifelse(d$AGEP_A >= 65 & d$AGEP_A <= 74,"65-74",
                   ifelse(d$AGEP_A >= 75 & d$AGEP_A <= 84,"75-84",
                   ifelse(d$AGEP_A == 85,"85+",NA_character_))),
                   levels=c("65-74","75-84","85+"))
d$sex <- fmap(d$SEX_A,c('1'='Male','2'='Female'),c('Male','Female'))
d$raceeth <- fmap(d$HISPALLP_A,
  c('1'='Hispanic','2'='Non-Hispanic White','3'='Non-Hispanic Black',
    '4'='Non-Hispanic Asian','5'='Other','6'='Other','7'='Other'),
  c('Non-Hispanic White','Non-Hispanic Black','Hispanic','Non-Hispanic Asian','Other'))
d$education <- fmap(d$EDUCP_A,
  c('1'='<High school/GED','2'='<High school/GED','3'='<High school/GED',
    '4'='High school graduate', '5'='Some college/associate',
    '6'='Some college/associate','7'='Some college/associate',
    '8'='Bachelor or higher','9'='Bachelor or higher','10'='Bachelor or higher'),
  c('Bachelor or higher','<High school/GED','High school graduate','Some college/associate'))
d$marital <- fmap(d$MARSTAT_A,
  c('1'='Married/partnered','2'='Married/partnered',
    '3'='Married/partnered','8'='Married/partnered',
    '4'='Not married/partnered','5'='Not married/partnered',
    '6'='Not married/partnered','7'='Not married/partnered'),
  c('Married/partnered','Not married/partnered'))
d$bmi <- fmap(d$BMICAT_A,
  c('1'='Underweight','2'='Healthy weight','3'='Overweight','4'='Obesity'),
  c('Healthy weight','Underweight','Overweight','Obesity'))
d$smoking <- fmap(d$SMKCIGST_A,
  c('1'='Current','2'='Current','3'='Former','4'='Never'),
  c('Never','Former','Current'))
d$phq2 <- fmap(d$PHQ2SCREEN_A,
  c('1'='Positive','2'='Negative'),c('Negative','Positive'))
d$vision_imp <- fmap(d$VISIONDF_A,
  c('1'='No','2'='Yes','3'='Yes','4'='Yes'),c('No','Yes'))
d$hearing_imp <- fmap(d$HEARINGDF_A,
  c('1'='No','2'='Yes','3'='Yes','4'='Yes'),c('No','Yes'))
for (pair in list(c('HYPEV_A','hypertension'), c('DIBEV_A','diabetes'),
                  c('COPDEV_A','copd'), c('CHDEV_A','chd'),
                  c('STREV_A','stroke'),c('BALDIZZ_A','balance'))) {
  d[[pair[2]]] <- fmap(d[[pair[1]]],c('1'='Yes','2'='No'),c('No','Yes'))
}
d$mobility <- fmap(d$DIFF_A,
  c('1'='No difficulty','2'='Any difficulty',
    '3'='Any difficulty','4'='Any difficulty'),
  c('No difficulty','Any difficulty'))

# 3. CONSTRUCT EACH MODEL'S EXACT COMPLETE-CASE SAMPLE ------------------------
vars1 <- c("recurrent", "pain_cat")
vars2 <- c(vars1,"agegrp","sex","raceeth","education","marital")
vars3 <- c(vars2,"bmi","smoking","phq2","vision_imp","hearing_imp",
           "hypertension","diabetes","copd","chd","stroke")
vars4 <- c(vars3,"balance","mobility")
for (k in 1:4) {
  vv <- get(paste0("vars",k))
  d[[paste0("analysis_m",k)]] <- as.integer(
    d$primary == 1 & complete.cases(d[,vv,drop=FALSE]))
}
d$complete4 <- d$analysis_m4

n_all <- nrow(d)
n_primary <- sum(d$primary)
n_falls <- sum(d$recurrent[d$primary == 1])
n_models <- sapply(1:4, function(k) sum(d[[paste0("analysis_m",k)]]))
cat(sprintf("Full dataset = %d, eligible = %d, recurrent events = %d\n",
            n_all, n_primary, n_falls))
cat("Model 1-4 sample sizes:", paste(n_models,collapse=", "),"\n")
if (n_all != 29522 || n_primary != 4473 || n_falls != 709 || n_models[4] != 4296) {
  stop("Sample verification failed: inspect recoding and NHIS file version first.")
}

# 4. NHIS DESIGN: DEFINE ON ALL 29,522 RECORDS, THEN DOMAIN SUBSET ------------
full_design <- svydesign(ids = ~PPSU, strata = ~PSTRAT, weights = ~WTFA_A,
                         data = d, nest = TRUE)
nominal_df <- degf(full_design) # bookkeeping = PSU count - stratum count
n_strata <- length(unique(d$PSTRAT))
n_psu <- nrow(unique(d[,c("PSTRAT","PPSU")]))
cat(sprintf("Full public-use design: %d strata, %d stratum-PSUs, degf=%d\n",
            n_strata,n_psu,nominal_df))
stopifnot(n_strata == 52, n_psu == 662, nominal_df == 610)

# Subsetting DESIGN instead of raw data preserves the original PSU/strata
# structure needed for domain-specific Taylor linearization.
design1 <- subset(full_design, analysis_m1 == 1)
design2 <- subset(full_design, analysis_m2 == 1)
design3 <- subset(full_design, analysis_m3 == 1)
design4 <- subset(full_design, analysis_m4 == 1)
designs <- list(design1,design2,design3,design4)
primary_design <- subset(full_design, primary == 1)

# 5. FORMULAS AND INDEPENDENT FITS --------------------------------------------
f1 <- recurrent ~ pain_cat
f2 <- update(f1, . ~ . + agegrp + sex + raceeth + education + marital)
f3 <- update(f2, . ~ . + bmi + smoking + phq2 + vision_imp + hearing_imp +
                          hypertension + diabetes + copd + chd + stroke)
f4 <- update(f3, . ~ . + balance + mobility)
formulas <- list(f1,f2,f3,f4)

# R's survey GLM estimates Taylor-linearized, design-adjusted SEs.
# quasi-binomial prevents meaningless integer-success warnings from weights.
run_fit <- function(formula, des) {
  svyglm(formula, design=des, family=quasibinomial(link="logit"),
         rescale=TRUE, na.action=na.omit)
}
fits <- lapply(1:4, function(k) run_fit(formulas[[k]], designs[[k]]))

# Extract coefficients and explicitly recompute inference under multiple
# reference distributions. df=120 is a sensitivity choice, NOT an NHIS rule.
extract_pain <- function(fit, model_name, sample_name, n, events, df=610) {
  b <- coef(fit); v <- vcov(fit); s <- sqrt(diag(v))
  targets <- c("pain_cat4-6","pain_cat7-10")
  if (!all(targets %in% names(b))) {
    stop("Pain coefficient names changed: ",paste(names(b),collapse=", "))
  }
  r_df <- tryCatch(as.numeric(summary(fit)$df[2]), error=function(e) NA_real_)
  if (length(r_df) != 1 || !is.finite(r_df) || r_df <= 0) r_df <- NA_real_
  result <- list()
  for (term in targets) {
    beta <- unname(b[term]); se <- unname(s[term]); z <- beta/se
    for (rule in c("t610_nominal","t120_sensitivity","normal","R_default_residual_df")) {
      dof <- switch(rule, t610_nominal=df, t120_sensitivity=120,
                    normal=Inf, R_default_residual_df=r_df)
      if (is.na(dof)) next
      cutoff <- if (is.infinite(dof)) qnorm(.975) else qt(.975,dof)
      pv <- if (is.infinite(dof)) 2*pnorm(abs(z),lower.tail=FALSE) else
              2*pt(abs(z), df=dof, lower.tail=FALSE)
      result[[length(result)+1]] <- data.frame(
        sample=sample_name, model=model_name, contrast=sub("pain_cat","",term),
        term=term, N=n, events=events, beta=beta, SE=se,
        OR=exp(beta), CI_low=exp(beta-cutoff*se),
        CI_high=exp(beta+cutoff*se), p=pv, inference=rule,
        reference_df=dof, stringsAsFactors=FALSE)
    }
  }
  do.call(rbind,result)
}

n_each <- sapply(1:4, function(k) sum(d[[paste0("analysis_m",k)]]))
e_each <- sapply(1:4, function(k) sum(d$recurrent[d[[paste0("analysis_m",k)]]==1]))
models_table <- do.call(rbind,lapply(1:4,function(k)
  extract_pain(fits[[k]],paste0("Model ",k),"original_samples",
               n_each[k],e_each[k],nominal_df)))
write.csv(models_table,file.path(out_dir,"R_original_models_all_df_options.csv"),row.names=FALSE)

# 6. SAME 4,296 PEOPLE IN EACH MODEL ------------------------------------------
fixed_fits <- lapply(formulas,function(f) run_fit(f,design4))
fixed_n <- sum(d$complete4)
fixed_e <- sum(d$recurrent[d$complete4==1])
fixed_table <- do.call(rbind,lapply(1:4,function(k)
  extract_pain(fixed_fits[[k]],paste0("Model ",k),"fixed_4296",
               fixed_n,fixed_e,nominal_df)))
write.csv(fixed_table,file.path(out_dir,"R_fixed_4296_models_all_df_options.csv"),row.names=FALSE)

# 7. INDEPENDENT PREVALENCE CHECK ---------------------------------------------
prev_rows <- list()
for (g in c("0-3","4-6","7-10")) {
  domain <- subset(primary_design, pain_cat == g)
  mean_result <- svymean(~recurrent, domain, na.rm=TRUE)
  estimate <- as.numeric(coef(mean_result)[1]); se <- as.numeric(SE(mean_result)[1])
  critical <- qt(.975,nominal_df)
  selection <- d$primary==1 & !is.na(d$pain_cat) & d$pain_cat==g
  prev_rows[[length(prev_rows)+1]] <- data.frame(
    pain=g, N=sum(selection), events=sum(d$recurrent[selection]),
    weighted_pct=100*estimate, SE_pct=100*se,
    low_pct=100*max(0,estimate-critical*se),
    high_pct=100*min(1,estimate+critical*se))
}
prev_table <- do.call(rbind,prev_rows)
write.csv(prev_table,file.path(out_dir,"R_prevalence_by_pain.csv"),row.names=FALSE)

# 8. DIRECT COMPARISON TO PYTHON'S FULL-PRECISION MODEL 4 RESULTS --------------
# The original Python's Model 4 point estimates were independently reproduced
# before this R run. This comparison does not assume their SEs are correct.
python_model4 <- data.frame(
  contrast=c("4-6","7-10"),
  Python_OR=c(1.2735382805821989,1.7471452669600710),
  Python_CI_low=c(1.0003671184285972,1.3356104295841964),
  Python_CI_high=c(1.6213045413327720,2.2854842372047743),
  Python_p_t610=c(0.04965380359220462,0.00005107037055789104))
# Recover full-precision Python log-odds SE from CI widths using original t610.
python_model4$Python_SE <-
  (log(python_model4$Python_CI_high)-log(python_model4$Python_CI_low)) /
  (2*qt(.975,610))
r4 <- subset(models_table,model=="Model 4" & inference=="t610_nominal")
check <- merge(python_model4,r4[,c("contrast","N","OR","SE","CI_low","CI_high","p")],
               by="contrast",sort=FALSE)
check$absolute_OR_difference <- abs(check$OR-check$Python_OR)
check$relative_SE_difference_pct <- 100*(check$SE/check$Python_SE-1)
check$absolute_p_difference <- abs(check$p-check$Python_p_t610)
write.csv(check,file.path(out_dir,"R_vs_Python_Model4_comparison.csv"),row.names=FALSE)

# Optional: if original Python CSV files exist, compare ALL Model 1-4 betas/SEs.
py_dir <- file.path(script_dir, "nhis_jointpain_project")
py_all <- list()
for (k in 1:4) {
  fn <- file.path(py_dir,paste0("model_",k,"_full_coefficients.csv"))
  if (!file.exists(fn)) next
  py <- read.csv(fn,stringsAsFactors=FALSE)
  py <- py[grepl("pain_cat",py$term,fixed=TRUE) &
             grepl("[T.",py$term,fixed=TRUE),]
  if (!nrow(py)) next
  py$contrast <- ifelse(grepl("4-6",py$term,fixed=TRUE),"4-6","7-10")
  py$model <- paste0("Model ",k)
  py_all[[length(py_all)+1]] <- py[,c("model","contrast","beta","se","OR")]
}
if (length(py_all)) {
  py <- do.call(rbind,py_all)
  names(py) <- c("model","contrast","Python_beta","Python_SE","Python_OR")
  rr <- subset(models_table,inference=="t610_nominal")
  joint <- merge(py,rr[,c("model","contrast","beta","SE","OR")],
                 by=c("model","contrast"))
  joint$beta_difference <- joint$beta-joint$Python_beta
  joint$SE_difference <- joint$SE-joint$Python_SE
  write.csv(joint,file.path(out_dir,"R_vs_Python_all_models_if_available.csv"),
            row.names=FALSE)
}

# 9. HUMAN-READABLE CONSOLE REPORT --------------------------------------------
cat("\n===== 2023 NHIS INDEPENDENT R survey VALIDATION =====\n")
cat("Original model sample sizes:",paste(n_each,collapse=", "),"\n")
cat("Common sample size / recurrent events:",fixed_n,"/",fixed_e,"\n")
cat("Bookkeeping degf(full_design):",nominal_df,"\n")
cat("degf(primary_design):",degf(primary_design),"\n")
cat("\nOR, 95% CI and p-values calculated with the same nominal t610 as Python:\n")
a <- subset(models_table,inference=="t610_nominal")
print(a[,c("model","contrast","N","events","OR","SE","CI_low","CI_high","p")],
      row.names=FALSE,digits=7)
cat("\nCommon 4,296-person sample, using t610:\n")
b <- subset(fixed_table,inference=="t610_nominal")
print(b[,c("model","contrast","N","OR","CI_low","CI_high","p")],
      row.names=FALSE,digits=7)
cat("\nPrevalence estimates by pain category:\n")
print(prev_table,row.names=FALSE,digits=7)
cat("\nDirect R vs Python Model 4 comparison:\n")
print(check[,c("contrast","Python_OR","OR","absolute_OR_difference",
               "Python_SE","SE","relative_SE_difference_pct","Python_p_t610","p")],
      row.names=FALSE,digits=9)
cat("\nIMPORTANT: Statistical/software agreement does NOT establish that\n")
cat("610 is the effective NHIS inferential df; NCHS cautions otherwise.\n")
cat("Inspect t120 and R-default results before interpreting p near 0.05.\n")
cat("CSV results saved to:",normalizePath(out_dir),"\n")
capture.output(sessionInfo(),file=file.path(out_dir,"R_sessionInfo.txt"))
