Main analyses for "Resting heart rate and antisocial behaviour: 
A Mendelian randomisation study."
================
Lucy Karwatowska

# Description

This is an Markdown document detailing the analysis steps for the paper titled ["Resting heart rate and antisocial behaviour: A Mendelian randomisation study."](link to manuscript). 


# Contents
-   [Load packages](#load-packages)
-   [Load outcome data](#load-outcome-data)
-   [Process outcome data](#process-outcome-data)
-   [Load exposure data](#load-exposure-data)
-   [Process exposure data](#process-exposure-data)
-   [Quality control procedures](#quality-control-procedures)
-         * [A. Select significant SNPs](A.-select-significant-SNPs)
-         * [B. Only keep exposure SNPs which are also in the outcome dataset](B.-only-keep-exposure-SNPs-which-are-also-in-the-outcome-dataset)
        - [C. Clump the exposure data](C.-clump-the-exposure-data)
        - [D. Harmonise the summary statistics](D.-harmonise-the-summary-statistics)
-   [Univariate Mendelian randomisation analysis between RHR and ASB](#univariate-mendelian-randomisation-analysis-between-rhr-an-asb)
-   [Sensitivity analyses](#sensitivity-analyses)
        - [A. Heterogeneity statistics](#A.-heterogeneity-statistics)
        - [B. Horizontal pleiotropy](#B.-horizontal-pleiotropy)
        - [C. Leave-one-out-analyses](#C.-leave-one-out-analyses)
        - [D. MR Egger](#D.-mr-egger)
        - [E. MR Steiger](#E.-mr-steiger)
        - [F. MR-PRESSO](#F.-mr-presso)
        - [G. Contamination mixture](#G.-contamination-mixture)
-   [Rerun analyses](#rerun-analyses)
-   [Plots](#plots)
        - [A. Scatter plot](#A.-scatter-plot)
        - [B. Forest plot](#A.-forest-plot)
        - [C. Funnel plot](#A.-funnel-plot)
-   [Reports](#reports)


      
# Load packages

``` r
library(arrangements) # 1.1.9
library(data.table) # 1.14.8
library(devtools) # 2.4.5
library(dpylr) # 1.1.2
library(gmodels) # 2.18.1.1
library(here) # 1.0.1
library(MendelianRandomization) # 0.7.0
library(mr.raps) # 0.2
library(MRPRESSO) # 1.0
library(TwoSampleMR) # 0.5.7

```


# Load outcome data

``` r

# Import the summary statistics file into R
SumStats <- fread(here::here("Data", "ASB_MR_Tielbeek_2021.txt"), header = T, data.table = T) 

head(SumStats)
tail(SumStats)

# The ASB summary statistics containing the converted effect size can now be saved as a new file in our "data" folder
write.table(SumStats, "Data/ASB_MR_Tielbeek_2021.txt", col.names = T, row.names = F, quote = F, sep = "\t")  


```


# Process outcome data

``` r
## read outcome data
outcome_ASB_dat <- read_outcome_data(filename = here::here("Data", "ASB_MR_Tielbeek_2021.txt"),
                                     snps = NULL, # Keep 'snps = NULL' so that all SNPs are uploaded
                                     sep = "\t",  # file tab ("\t") or  space delimited ("")
                                     snp_col = "SNP",  # Insert the name of column that includes the SNP IDs
                                     beta_col = "BETA", # Insert the name of column with effect sizes
                                     se_col = "SE", # Insert the name of column with standard errors
                                     eaf_col = "MAF",  # effect allele frequency
                                     effect_allele_col = "A1", # Insert the name of column with the effect allele
                                     other_allele_col = "A2",  # Insert the name of column with the non-effect allele
                                     pval_col = "P",# Insert the name of column with p-value
                                     samplesize_col = "N") # sample size available for each SNP

# Note: A few columns are not present in the sum stat file so cannot be read directly but can be added

# Provide a name for the outcome
outcome_ASB_dat$outcome <- "ASB" 

head(outcome_ASB_dat)

# How many SNPs are included?
length(outcome_ASB_dat$SNP) # Number of SNPs included 6565696

# Check the number of rows with command line
system("wc -l Data/ASB_MR_Tielbeek_2021.txt") # one more than the R object as previously as counts column names as one row 
system("head Data/ASB_MR_Tielbeek_2021.txt")
system("tail Data/ASB_MR_Tielbeek_2021.txt")
```


# Load exposure data

``` r
# Read resting heart rate as exposure data

## Import the summary statistics file into R
SumStatsRestingHR <- fread(here::here("Data", "RHRnoSmoke_MR_Zhu_2019.txt"), header = T, data.table = T) 

head(SumStatsRestingHR)

# The summary statistics containing the converted effect size can now be exported and saved as a new file in our "data" folder
write.table(SumStatsRestingHR, "Data/RHRnoSmoke_MR_Zhu_2019.txt", col.names = T, row.names = F, quote = F, sep = "\t") 
```


# Process exposure data

``` r
# read data using dedicated function of the TwoSampleMR package (pendent to the read_outcome_data() function)
exposure_HR_dat <- read_exposure_data(filename = "Data/RHRnoSmoke_MR_Zhu_2019.txt", 
                                      sep = "\t",
                                      snp_col = "SNP",           
                                      beta_col = "BETA",         
                                      se_col = "SE",       
                                      effect_allele_col = "A1",  
                                      other_allele_col = "A2",   
                                      pval_col = "P",  
                                      eaf_col = "MAF") 

# Provide a name for the outcome  
exposure_HR_dat$exposure <- "RestingHeartRate"     

# Add the sample size manually from the publication, as not available in the summary statistics
exposure_HR_dat$samplesize.exposure <- 458835

# How many SNPs?
length(exposure_HR_dat$SNP) # N = 9705784 SNPs
```


# Quality control procedures

## A: Select significant SNPs

``` r
# Select genome-wide significant SNPs
exposure_HR_dat5e8 <- subset(exposure_HR_dat, pval.exposure < 5e-8)


# Lets have a look at the number of SNPs that we would include:
length(exposure_HR_dat5e8$SNP) # 46766

```


## B: Only keep exposure SNPs which are also in the outcome dataset 

``` r
# Keep all SNPs from the exposure set which are present in the outcome set
exposure_HR_dat5e8_out <- exposure_HR_dat5e8[exposure_HR_dat5e8$SNP %in% outcome_ASB_dat$SNP, ]

# Lets have a look at the number of SNPs we excluded
length(exposure_HR_dat5e8_out$SNP) - length(exposure_HR_dat5e8$SNP) # N excluded -18087

```


## C: Clump the exposure data

``` r
# using default (recommended) parameters
exposure_HR_dat5e8_clumped <- clump_data(exposure_HR_dat5e8_out) # Using default parameters

# Lets have a look at the number of SNPs we excluded
length(exposure_HR_dat5e8_clumped$SNP) # N included = 300
length(exposure_HR_dat5e8_clumped$SNP) - length(exposure_HR_dat5e8_out$SNP) # N excluded - 28379


# F statistic
# ((n-k-1)/k) * (R2/ 1-R2)
# ((458835-297-1)/297) * (0.01 / (1 - 0.01))
# ((458835-297-1)/297) * (0.02 / (1 - 0.02)) 
# ((458835-297-1)/297) * (0.05 / (1 - 0.05)) 
```


## D: Harmonise the summary statistics

``` r
####### Harmonize
DataMR_HRtoASB <- harmonise_data(exposure_dat = exposure_HR_dat5e8_clumped,
                                 outcome_dat = outcome_ASB_dat, 
                                 action = 2)

table(DataMR_HRtoASB$mr_keep) 
# This columns tells you which SNPs will be kept for the main analysis: All rows that are set to TRUE will be included in the MR analysis. 

attr(DataMR_HRtoASB, "log") # Detailed summary of what was done and reasons for excluding SNPs 

# F statistics

DataMR_HRtoASB$fstatistic <- ((DataMR_HRtoASB$beta.exposure^2) / (DataMR_HRtoASB$se.exposure^2))

nrow(DataMR_HRtoASB[DataMR_HRtoASB$fstatistic>10, ]) # 300 (so all!)
mean(DataMR_HRtoASB$fstatistic) # 81.34597
summary(DataMR_HRtoASB$fstatistic)
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 27.30   36.25   47.19   81.35   73.41 1185.15 
```


# Univariate Mendelian randomisation analysis between RHR and ASB

``` r

# Run as subset of estimators
(MR_HRtoASB <- mr(DataMR_HRtoASB, 
                  method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression")))

knitr::kable(as.data.frame(MR_HRtoASB), "markdown")


# mr_raps

mr_frame1 <- as.data.frame(DataMR_HRtoASB[DataMR_HRtoASB$mr_keep == TRUE, # keep only SNPS after harmoization
                                           c('beta.exposure', 'se.exposure', 'beta.outcome', 'se.outcome')])
names(mr_frame1) <- c("RHR","RHR_se","ASB","ASB_se")

# Perform MR-RAPS 
(mr_raps_result <- mr.raps.overdispersed.robust(b_exp = mr_frame1$RHR, b_out = mr_frame1$ASB, # define effects of SNPs on exposure and outcome
                             se_exp = mr_frame1$RHR_se, se_out = mr_frame1$ASB_se, # define standard errors of exposure and outcome
                             loss.function = "huber", # choose loss function from c("huber", "tukey")
                             k = 1.345, # Threshold parameter in the Huber and Tukey loss functions. (default)
                             initialization = c("l2"), # Method to initialize the robust estimator, c("l2", "mode")
                             suppress.warning = FALSE, 
                             niter = 20, # Maximum number of interations to solve the estimating equations. (default)
                             tol = .Machine$double.eps^0.5)) # Numerical precision. (default)

knitr::kable(as.data.frame(mr_raps_result), "markdown")



```


``` r
write.csv(MR_HRtoASB, file = "output/Main/MR_restingHRtoASB_unfiltered.csv") # Save the results 

```


# Sensitivity analyses

## A. Heterogeneity statistics

``` r

# Cochran Q statistics for heterogeneity using default methods (mr_egger and mr_ivw)
(Results_HRtoASB_heterogeneity <- mr_heterogeneity(DataMR_HRtoASB))
                      
knitr::kable(as.data.frame(Results_HRtoASB_heterogeneity), "markdown")

# results indicate significant heterogeneity of SNPs -> potentially due to pleiotropy

```


## B. Horizontal pleiotropy  

``` r
# Pleiotropy test
(Results_HRtoASB_pleiotropy <- mr_pleiotropy_test(DataMR_HRtoASB))

knitr::kable(as.data.frame(Results_HRtoASB_pleiotropy), "markdown")

# no significant directional pleiotropy detected  p = 0.985006
```


## C. Leave-one-out analysis

``` r

# Leave-one-out analysis
(DataMR_HRtoASB_leaveOut_ivw <- mr_leaveoneout(DataMR_HRtoASB)) # Default: ivw method
DataMR_HRtoASB_leaveOut_egger <- mr_leaveoneout(DataMR_HRtoASB, method = mr_egger_regression)


knitr::kable(as.data.frame(DataMR_HRtoASB_leaveOut_ivw), "markdown")
knitr::kable(as.data.frame(DataMR_HRtoASB_leaveOut_egger), "markdown")
```


``` r
# Check max and min of the beta and the pvalue to see if any snp has changed them dramatically
out1 <- apply(DataMR_HRtoASB_leaveOut_ivw[ , c("b","p")], 2, 
              function(x) {cbind(min(x), max(x))})

knitr::kable(as.data.frame(out1), "markdown")

# See that estimates vary 
out2 <- apply(DataMR_HRtoASB_leaveOut_egger[ , c("b","p")], 2, 
              function(x) {cbind(min(x), max(x))})

knitr::kable(as.data.frame(out2), "markdown")


```

## D. MR Egger

``` r
# Calculate baselines estimates on this new dataset with weaker instruments
# Preliminary step: take only variants with smallest absolute effect sizes, if not (we have I2 >99%) and no correction needed

# First, select variants with the weaker betas
DataMR_HRtoASBk <- DataMR_HRtoASB[abs(DataMR_HRtoASB$beta.exposure) < quantile(abs(DataMR_HRtoASB$beta.exposure), 0.20), ] 

# Second, take out instruments excluded in  TwoSampleMR
DataMR_HRtoASBk_keep <- DataMR_HRtoASBk[DataMR_HRtoASBk$mr_keep == "TRUE", ] 

length(DataMR_HRtoASBk_keep$SNP) # 54 variants selected.
# Note that a very small number for MR Egger (e.g. 8) yield completely unstable estimates.

# Calculate baselines estimates on this new dataset with weaker instruments
(MR_HRtoASBnew <- mr(DataMR_HRtoASBk_keep, 
                     method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression")))

knitr::kable(as.data.frame(MR_HRtoASBnew), "markdown")

# no large differences when using only weak instruments compared to the original results (still no causal effect)


```


``` r
#Function to calculate I2
Isq <- function(y,s){
  k = length(y)
  w = 1/s^2; sum.w = sum(w)
  mu.hat = sum(y*w)/sum.w
  Q = sum(w*(y - mu.hat)^2)
  Isq = (Q - (k - 1))/Q
  Isq = max(0, Isq)
  return(Isq)
}


```


``` r
# 1. Calculate I2 (I squared) - include only weak instruments
(isq_w <- Isq(DataMR_HRtoASBk_keep$beta.exposure / DataMR_HRtoASBk_keep$se.outcome,
                 DataMR_HRtoASBk_keep$se.exposure / DataMR_HRtoASBk_keep$se.outcome)) 

# Despite weaker variants, I2 is high so really no problem of weak instrument with resting heart rate
# 0.9716064
```


``` r
# 2. Crude correction BMR/I2 (NOT required to correct because Isquare very high):
mr(DataMR_HRtoASB, method_list = c("mr_egger_regression"))
mr(DataMR_HRtoASB, method_list = c("mr_egger_regression"))$b/isq

# inference stays the same, even when correcting
```

## E. MR Steiger


``` r
## MR Steiger
#     - TwoSampleMR relies on package psych, r.test(n = mean(n_exp), n2 = mean(n_out), r12 = r_exp, r34 = r_out)
#     -  n and n2 are the average N for exposure and outcome, r12 the r of the correlation with exposure and r34 with outcome (the 12 and 34 are there two show that variables are not the same)

# 2) Run MR-Steiger
(DataMR_HRtoASB_Steiger <- directionality_test(DataMR_HRtoASB))

knitr::kable(as.data.frame(DataMR_HRtoASB_Steiger), "markdown")

DataMR_HRtoASB_Steiger$snp_r2.exposure / DataMR_HRtoASB_Steiger$snp_r2.outcome # effect on the exposure 7.485474 times bigger

# correct direction as expected (R2 for the outcome is super small!)
```

``` r
# We therefore apply the test on one SNP
DataMR_HRtoASB_Steiger <- directionality_test(DataMR_HRtoASB[1, ])

# As can be seen the r2 for the exposure is much bigger than for the outcome (although both of them are very small as this is just one SNP)
DataMR_HRtoASB_Steiger$snp_r2.exposure / DataMR_HRtoASB_Steiger$snp_r2.outcome # effect on the exposure 2.3 times bigger
```


``` r

# Loop over all SNPs
SteigerSNPs = matrix(nrow = dim(DataMR_HRtoASB)[1], ncol = 8) 

for (i in 1:dim(DataMR_HRtoASB)[1]){
  temp  <- directionality_test(DataMR_HRtoASB[i, ])
  SteigerSNPs[i, ] = as.matrix(temp)
}

SteigerSNPs <- as.data.frame(SteigerSNPs)
names(SteigerSNPs) <- names(DataMR_HRtoASB_Steiger)

knitr::kable(head(SteigerSNPs), "markdown")
```


``` r
table(SteigerSNPs$correct_causal_direction) # most SNPS valid instruments, 14 FALSE

wrongDirect <- DataMR_HRtoASB[which(SteigerSNPs$correct_causal_direction == F), "SNP"]
```


``` r
DataMR_HRtoASB$r.exposure <- get_r_from_pn(n = DataMR_HRtoASB$samplesize.exposure,
                                           p = DataMR_HRtoASB$pval.exposure)

DataMR_HRtoASB$r.outcome <- get_r_from_pn(n = DataMR_HRtoASB$samplesize.outcome,
                                          p = DataMR_HRtoASB$pval.outcome)

# and recalculate Steiger
SteigerSNPs1 = matrix(nrow = dim(DataMR_HRtoASB)[1], ncol = 8) 

for (i in 1:dim(DataMR_HRtoASB)[1]){
  temp  <- directionality_test(DataMR_HRtoASB[i, ])
  SteigerSNPs1[i, ] = as.matrix(temp)
}

SteigerSNPs1 = as.data.frame(SteigerSNPs1)

names(SteigerSNPs1) = names(DataMR_HRtoASB_Steiger)  
table(SteigerSNPs1$correct_causal_direction) # most SNPS valid instruments

head(SteigerSNPs1)

# As expected, these are identical as we are just using the internal function used by directionality_test() to calculate r
identical(SteigerSNPs$steiger_pval, SteigerSNPs1$steiger_pval)
```


### Rerun the univariate MR analyses between RHR and ASB using only the appropriate instruments

- In other words, without the alleles identified in the MR steiger filtering

``` r
# integrate the correct_causal_direction column into DataMR
DataMR_HRtoASB$correct_causal_direction <- SteigerSNPs1$correct_causal_direction 

# create DataMR only with those SNPs
DataMR_HRtoASB_Steigerfiltered <- DataMR_HRtoASB[DataMR_HRtoASB$correct_causal_direction == "TRUE", ] 

# check that number of rows matches
dim(DataMR_HRtoASB_Steigerfiltered)

# Redo MR using only valid instruments (SNPs)
(MR_HRtoASB_Steigerfiltered <- mr(DataMR_HRtoASB_Steigerfiltered,
                     method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression")))

knitr::kable(MR_HRtoASB_Steigerfiltered, "markdown")

write.csv(MR_HRtoASB_Steigerfiltered, file = "output/Main/MR_restingHRtoASB_SteigerFiltered.csv") # Saves the results 

# no large deviations from the original results (still no causal effect)


# mr_raps

mr_frame2 <- as.data.frame(DataMR_HRtoASB_Steigerfiltered[DataMR_HRtoASB_Steigerfiltered$mr_keep == TRUE, # keep only SNPS after harmonization
                                           c('beta.exposure', 'se.exposure', 'beta.outcome', 'se.outcome')])
names(mr_frame2) <- c("RHR","RHR_se","ASB","ASB_se")

# Perform MR-RAPS 
(mr_raps_result2 <- mr.raps.overdispersed.robust(b_exp = mr_frame2$RHR, b_out = mr_frame2$ASB, # define effects of SNPs on exposure and outcome
                             se_exp = mr_frame2$RHR_se, se_out = mr_frame2$ASB_se, # define standard errors of exposure and outcome
                             loss.function = "huber", # choose loss function from c("huber", "tukey")
                             k = 1.345, # Threshold parameter in the Huber and Tukey loss functions. (default)
                             initialization = c("l2"), # Method to initialize the robust estimator, c("l2", "mode")
                             suppress.warning = FALSE, 
                             niter = 20, # Maximum number of interations to solve the estimating equations. (default)
                             tol = .Machine$double.eps^0.5)) # Numerical precision. (default)

knitr::kable(as.data.frame(mr_raps_result2), "markdown")


# CIs for different estimates
(CI_lower_ivw <- MR_HRtoASB_Steigerfiltered$b[1] - 1.96*MR_HRtoASB_Steigerfiltered$se[1])
(CI_upper_ivw <- MR_HRtoASB_Steigerfiltered$b[1] + 1.96*MR_HRtoASB_Steigerfiltered$se[1])

(CI_lower_weightedM <- MR_HRtoASB_Steigerfiltered$b[2] - 1.96*MR_HRtoASB_Steigerfiltered$se[2])
(CI_upper_weightedM <- MR_HRtoASB_Steigerfiltered$b[2] + 1.96*MR_HRtoASB_Steigerfiltered$se[2])

(CI_lower_egg <- MR_HRtoASB_Steigerfiltered$b[3] - 1.96*MR_HRtoASB_Steigerfiltered$se[3])
(CI_upper_egg <- MR_HRtoASB_Steigerfiltered$b[3] + 1.96*MR_HRtoASB_Steigerfiltered$se[3])

(CI_lower_raps <- mr_raps_result2$beta.hat - 1.96*mr_raps_result2$beta.se)
(CI_upper_raps <- mr_raps_result2$beta.hat + 1.96*mr_raps_result2$beta.se)

# also repeat heterogeneity test using only valid instrumentsv

# Cochran Q statistics for heterogeneity using default methods (mr_egger and mr_ivw)
(Results_HRtoASB_heterogeneity_filtered <- mr_heterogeneity(DataMR_HRtoASB_Steigerfiltered))
                      
knitr::kable(as.data.frame(Results_HRtoASB_heterogeneity_filtered), "markdown")

# YES, no more significant heterogeneity!


# Pleiotropy test
(Results_HRtoASBsteiger_pleiotropy <- mr_pleiotropy_test(DataMR_HRtoASB_Steigerfiltered))

knitr::kable(as.data.frame(Results_HRtoASBsteiger_pleiotropy), "markdown")


# F statistics

DataMR_HRtoASB_Steigerfiltered$fstatistic <- ((DataMR_HRtoASB_Steigerfiltered$beta.exposure^2) / (DataMR_HRtoASB_Steigerfiltered$se.exposure^2))

DataMR_HRtoASB_Steigerfiltered <- DataMR_HRtoASB_Steigerfiltered[DataMR_HRtoASB_Steigerfiltered$mr_keep == TRUE, ]

nrow(DataMR_HRtoASB_Steigerfiltered[DataMR_HRtoASB_Steigerfiltered$fstatistic>10, ]) # 278 (so all!)
mean(DataMR_HRtoASB_Steigerfiltered$fstatistic) # 84.20383
summary(DataMR_HRtoASB_Steigerfiltered$fstatistic)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   27.30   36.53   48.34   84.20   76.69 1185.15 
```


## F. MR-PRESSO 


``` r

# requires data to be in a data frame
mr_frame <- as.data.frame(DataMR_HRtoASB_Steigerfiltered[DataMR_HRtoASB_Steigerfiltered$mr_keep == TRUE, # keep only SNPS after harmonization
                                         c('beta.exposure', 'se.exposure', 'beta.outcome', 'se.outcome')])

names(mr_frame) <- c("HR","HR_se","ASB","ASB_se")

# Perform MR-PRESSO
out_presso <- mr_presso(BetaOutcome = "ASB",   # define effects of SNPs on exposure and outcome
                        BetaExposure = "HR",   
                        SdOutcome = "ASB_se",  # define standard errors of exposure and outcome
                        SdExposure = "HR_se",  
                        OUTLIERtest = TRUE,    # correction of pleiotropy via outlier removal
                        DISTORTIONtest = TRUE, # testing of sign. distortion in causal estimate before & after MR-PRESSO correction
                        data = mr_frame, 
                        NbDistribution = 1000, # number of elements to simulate to form null distribution to compute empirical Ps
                        SignifThreshold = 0.05,
                        seed = 123)            # seed for reproducibility of estimates

out_presso



# no pleiotropy


# check the main results - raw vs corrected estimates 
res1 <- out_presso$`Main MR results`
knitr::kable(res1)

# confirming previous results of no causal effect 
# no outlier correction because no outliers were present

out_presso$`Main MR results`$`Causal Estimate`[1]
(CI_lower <- out_presso$`Main MR results`$`Causal Estimate`[1] - 1.96 * out_presso$`Main MR results`$Sd[1])
(CI_upper <- out_presso$`Main MR results`$`Causal Estimate`[1] + 1.96* out_presso$`Main MR results`$Sd[1])

# save output 
write.csv(res1, file = "output/Main/MRpresso_mainResults.csv", row.names = F, quote = F)
write.csv(out_presso$`MR-PRESSO results`, file = "output/Main/MRpresso_tests.csv", row.names = F, quote = F)


```


## G. Contamination mixture method 


``` r

mr_object2 <- mr_input(bx = DataMR_HRtoASB_Steigerfiltered$beta.exposure, 
                       bxse = DataMR_HRtoASB_Steigerfiltered$se.exposure,
                       by = DataMR_HRtoASB_Steigerfiltered$beta.outcome,
                       byse = DataMR_HRtoASB_Steigerfiltered$se.outcome)

# Perform Contamination 
(out_contam2 <- mr_conmix(mr_object2, 
                        psi = 0,       # psi = value of the standard deviation parameter (default = 0)
                        CIMin = NA, 
                        CIMax = NA, 
                        CIStep = 0.01, # The step size used in the search to find the confidence interval.
                        alpha = 0.05)) # significance level used when calculating the confidence intervals.

# N SNPs 278
out_contam2@Estimate # close to 0 -> -0.002325669
out_contam2@Psi # 0.05169542 (SD) 
se <- out_contam2@Psi/sqrt(out_contam2@SNPs)
out_contam2@Pvalue   #  0.4519773

(CI_lower <- out_contam2@Estimate - 1.96*se)
(CI_upper <- out_contam2@Estimate + 1.96*se)

```


# Plots

## A. Scatter plot


``` r
# Create scatter plot
MR_HRtoASBplots <- mr(DataMR_HRtoASB_Steigerfiltered, 
                                  method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))

MR_HRtoASB_scatter <- mr_scatter_plot(MR_HRtoASBplots, DataMR_HRtoASB_Steigerfiltered)
MR_HRtoASB_scatter[[1]]

# It is possible to save this plot using the ggsave() function, e.g. to save as a pdf or png
library(ggplot2)
ggplot2::ggsave(MR_HRtoASB_scatter[[1]], file = "output/Main/MR_restingHRtoASB_scatter_new.png", width = 7, height = 7)
```



## B. Forest plot

``` r

# Get the causal estimates for each SNP separately
MR_HRtoASB_single <- mr_singlesnp(DataMR_HRtoASB_Steigerfiltered, 
                                  all_method = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))
# Create forest plot
MR_HRtoASB_forest <- mr_forest_plot(MR_HRtoASB_single)
MR_HRtoASB_forest[[1]]

# save figure
ggplot2::ggsave(MR_HRtoASB_forest[[1]], file = "output/Main/MR_restingHRtoASB_forest_new.pdf", width = 15, height = 20)
```


## C. Funnel plot



``` r

# Create funnel plot using single SNPs
MR_HRtoASB_funnel <- mr_funnel_plot(MR_HRtoASB_single)
MR_HRtoASB_funnel[[1]]

# save figure
ggplot2::ggsave(MR_HRtoASB_funnel[[1]], file = "output/Main/MR_restingHRtoASB_funnel_new.png", width = 7, height = 7)
```


# Reports

``` r
# Report output
mr_report(DataMR_HRtoASB,  
          output_path = here::here("Output/Main"),
          study = "Mendelian Randomisation analysis: Resting heart rate and antisocial behaviour", 
          output_type = "html")


```
