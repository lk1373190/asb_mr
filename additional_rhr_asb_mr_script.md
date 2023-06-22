Additional analyses for "Resting heart rate and antisocial behaviour: 
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
        -  [A. Select significant SNPs](A.-select-significant-SNPs)
        -  [B. Only keep exposure SNPs which are also in the outcome dataset](B.-only-keep-exposure-SNPs-which-are-also-in-the-outcome-dataset)
        - [C. Clump the exposure data](C.-clump-the-exposure-data)
        - [D. Harmonise the summary statistics](D.-harmonise-the-summary-statistics)
-   [Univariate Mendelian randomisation analysis between RHR and ASB](#univariate-mendelian-randomisation-analysis-between-rhr-an-asb)
-   [Sensitivity analyses](#sensitivity-analyses)
        - [A. Heterogeneity statistics](#A.-heterogeneity-statistics)
        - [B. Horizontal pleiotropy](#B.-horizontal-pleiotropy)
        - [C. Leave-one-out-analyses](#C.-leave-one-out-analyses)
        - [D. MR Egger](#D.-mr-egger)
        - [E. MR Steiger](#E.-mr-steiger)
        - 


      
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

# Additional analyses

## Univariable MR analyses between HRV and ASB

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
```

### A. pvRSAHF


``` r
# Read HRV as exposure data

## Import the summary statistics file into R
SumStatsHRV <- fread(here::here("Data", 'HeartRateVarpvRSAHF_MR_Nolte_2017.txt') , header = T, data.table = T) 

head(SumStatsHRV)

# The summary statistics containing the converted effect size can now be exported and saved as a new file in our "data" folder
write.table(SumStatsHRV, "Data/HeartRateVarpvRSAHF_MR_Nolte_2017.txt", col.names = T, row.names = F, quote = F, sep = "\t") 

# read data using dedicated function of the TwoSampleMR package (pendent to the read_outcome_data() function)
exposure_HRV_pvRSAHF_dat <- read_exposure_data(filename = "Data/HeartRateVarpvRSAHF_MR_Nolte_2017.txt", 
                                      sep = "\t",
                                      snp_col = "SNP",           
                                      beta_col = "BETA",         
                                      se_col = "SE",       
                                      effect_allele_col = "A1",  
                                      other_allele_col = "A2",   
                                      pval_col = "P",  
                                      eaf_col = "MAF",
                                      samplesize_col = "N") 

# Provide a name for the outcome  
exposure_HRV_pvRSAHF_dat$exposure <- "HRV_pvRSAHF"     


# How many SNPs?
length(exposure_HRV_pvRSAHF_dat$SNP) # N = 2530479 SNPs
```


``` r
# Select SNPs for exposure

# Select genome-wide significant SNPs
exposure_HRV_pvRSAHF_dat5e8 <- subset(exposure_HRV_pvRSAHF_dat, pval.exposure < 5e-8)


# Select a more liberal threshold
exposure_HRV_pvRSAHF_dat5e5 <- subset(exposure_HRV_pvRSAHF_dat, pval.exposure < 5e-5)

# Lets have a look at the number of SNPs that we would include:
length(exposure_HRV_pvRSAHF_dat5e8$SNP) # 230
length(exposure_HRV_pvRSAHF_dat5e5$SNP) # 1080

```



``` r}
# Keep all SNPs from the exposure set which are present in the outcome set
exposure_HRV_pvRSAHF_dat5e8_out <- exposure_HRV_pvRSAHF_dat5e8[exposure_HRV_pvRSAHF_dat5e8$SNP %in% outcome_ASB_dat$SNP, ]

exposure_HRV_pvRSAHF_dat5e5_out <- exposure_HRV_pvRSAHF_dat5e5[exposure_HRV_pvRSAHF_dat5e5$SNP %in% outcome_ASB_dat$SNP, ]


# Lets have a look at the number of SNPs we excluded
length(exposure_HRV_pvRSAHF_dat5e8_out$SNP) - length(exposure_HRV_pvRSAHF_dat5e5_out$SNP) # N excluded 242

```



``` r
# Clump the data

# using default (recommended) parameters
exposure_HRV_pvRSAHF_dat5e8_clumped <- clump_data(exposure_HRV_pvRSAHF_dat5e8_out)

# Lets have a look at the number of SNPs we excluded
length(exposure_HRV_pvRSAHF_dat5e8_clumped$SNP) # N included = 4
length(exposure_HRV_pvRSAHF_dat5e8_clumped$SNP) - length(exposure_HRV_pvRSAHF_dat5e8_out$SNP) # N excluded - 55


# second threshold

exposure_HRV_pvRSAHF_dat5e5_clumped <- clump_data(exposure_HRV_pvRSAHF_dat5e5_out)

# Lets have a look at the number of SNPs we excluded
length(exposure_HRV_pvRSAHF_dat5e5_clumped$SNP) # N included = 44

```


``` r
# Harmonise summary statistics
DataMR_HRVpvRSAHFtoASB <- harmonise_data(exposure_dat = exposure_HRV_pvRSAHF_dat5e8_clumped,
                                 outcome_dat = outcome_ASB_dat, 
                                 action = 2)

table(DataMR_HRVpvRSAHFtoASB$mr_keep) #TRUE 4

# This columns tells you which SNPs will be kept for the main analysis: All rows that are set to TRUE will be included in the MR analysis. 

attr(DataMR_HRVpvRSAHFtoASB, "log") # Detailed summary of what was done and reasons for excluding SNPs 


# second threshold
DataMR_HRVpvRSAHFtoASB_lib <- harmonise_data(exposure_dat = exposure_HRV_pvRSAHF_dat5e5_clumped,
                                 outcome_dat = outcome_ASB_dat, 
                                 action = 2)

table(DataMR_HRVpvRSAHFtoASB_lib$mr_keep) # FALSE 1 TRUE 43


# F statistics

DataMR_HRVpvRSAHFtoASB$fstatistic <- ((DataMR_HRVpvRSAHFtoASB$beta.exposure^2) / (DataMR_HRVpvRSAHFtoASB$se.exposure^2))

nrow(DataMR_HRVpvRSAHFtoASB[DataMR_HRVpvRSAHFtoASB$fstatistic>10, ]) # 4 (so all!)
mean(DataMR_HRVpvRSAHFtoASB$fstatistic) # 86.01942
summary(DataMR_HRVpvRSAHFtoASB$fstatistic)
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 33.42   44.98   72.82   86.02  113.86  165.02 


# liberal threshold

DataMR_HRVpvRSAHFtoASB_lib <- DataMR_HRVpvRSAHFtoASB_lib[DataMR_HRVpvRSAHFtoASB_lib$mr_keep == TRUE, ]
DataMR_HRVpvRSAHFtoASB_lib$fstatistic <- ((DataMR_HRVpvRSAHFtoASB_lib$beta.exposure^2) / (DataMR_HRVpvRSAHFtoASB_lib$se.exposure^2))

nrow(DataMR_HRVpvRSAHFtoASB_lib[DataMR_HRVpvRSAHFtoASB_lib$fstatistic>10, ]) # 43 (so all!)
mean(DataMR_HRVpvRSAHFtoASB_lib$fstatistic) # 25.00546
summary(DataMR_HRVpvRSAHFtoASB_lib$fstatistic)
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #  16.49   17.08   18.50   25.01   20.62  165.02  
```


``` r

# Perform standard MR analyses
(MR_HRVpvRSAHFtoASB <- mr(DataMR_HRVpvRSAHFtoASB, 
                  method_list = c("mr_weighted_median","mr_egger_regression",
                                  "mr_ivw")))

knitr::kable(as.data.frame(MR_HRVpvRSAHFtoASB), "markdown")


# more liberal threshold
(MR_HRVpvRSAHFtoASB_lib <- mr(DataMR_HRVpvRSAHFtoASB_lib, 
                  method_list = c("mr_weighted_median", "mr_egger_regression",
                                  "mr_ivw")))

```


If you want to export your results in an excel file, you can run the following in R:
``` r
write.csv(MR_HRVpvRSAHFtoASB, file = "output/Additional/MR_HRV_pvRSAHFtoASB_unfiltered.csv") # Saves the results 
write.csv(MR_HRVpvRSAHFtoASB_lib, file = "output/Additional/MR_HRV_pvRSAHFtoASBlib_unfiltered.csv") # Saves the results 

```


### B. SDNN


``` r
# Read HRV as exposure data

## Import the summary statistics file into R
SumStatsHRV <- fread(here::here("Data", 'HeartRateVarSDNN_MR_Nolte_2017.txt') , header = T, data.table = T) 

head(SumStatsHRV)

# The summary statistics containing the converted effect size can now be exported and saved as a new file in our "data" folder
write.table(SumStatsHRV, "Data/HeartRateVarSDNN_MR_Nolte_2017.txt", col.names = T, row.names = F, quote = F, sep = "\t") 

# read data using dedicated function of the TwoSampleMR package (pendent to the read_outcome_data() function)
exposure_HRV_SDNN_dat <- read_exposure_data(filename = "Data/HeartRateVarSDNN_MR_Nolte_2017.txt", 
                                      sep = "\t",
                                      snp_col = "SNP",           
                                      beta_col = "BETA",         
                                      se_col = "SE",       
                                      effect_allele_col = "A1",  
                                      other_allele_col = "A2",   
                                      pval_col = "P",  
                                      eaf_col = "MAF",
                                      samplesize_col = "N") 

# Provide a name for the outcome  
exposure_HRV_SDNN_dat$exposure <- "HRV_SDNN"     


# How many SNPs?
length(exposure_HRV_SDNN_dat$SNP) # N = 2553521 SNPs
```


``` r
# Select SNPs for exposure

# Select genome-wide significant SNPs
exposure_HRV_SDNN_dat5e8 <- subset(exposure_HRV_SDNN_dat, pval.exposure < 5e-8)


# Select a more liberal threshold
exposure_HRV_SDNN_dat5e5 <- subset(exposure_HRV_SDNN_dat, pval.exposure < 5e-5)

# Lets have a look at the number of SNPs that we would include:
length(exposure_HRV_SDNN_dat5e8$SNP) # 103
length(exposure_HRV_SDNN_dat5e5$SNP) # 699

```


``` r
# Keep all SNPs from the exposure set which are present in the outcome set
exposure_HRV_SDNN_dat5e8_out <- exposure_HRV_SDNN_dat5e8[exposure_HRV_SDNN_dat5e8$SNP %in% outcome_ASB_dat$SNP, ]

exposure_HRV_SDNN_dat5e5_out <- exposure_HRV_SDNN_dat5e5[exposure_HRV_SDNN_dat5e5$SNP %in% outcome_ASB_dat$SNP, ]


# Lets have a look at the number of SNPs we excluded
length(exposure_HRV_SDNN_dat5e8_out$SNP) - length(exposure_HRV_SDNN_dat5e5_out$SNP) # N excluded 449

```


``` r
# Clump the data

# using default (recommended) parameters
exposure_HRV_SDNN_dat5e8_clumped <- clump_data(exposure_HRV_SDNN_dat5e8_out)

# Lets have a look at the number of SNPs we excluded
length(exposure_HRV_SDNN_dat5e8_clumped$SNP) # N included = 5
length(exposure_HRV_SDNN_dat5e8_clumped$SNP) - length(exposure_HRV_SDNN_dat5e8_out$SNP) # N excluded - 65


# second threshold

exposure_HRV_SDNN_dat5e5_clumped <- clump_data(exposure_HRV_SDNN_dat5e5_out)

# Lets have a look at the number of SNPs we excluded
length(exposure_HRV_SDNN_dat5e5_clumped$SNP) # N included = 54

```


``` r
# Harmonise summary statistics
DataMR_HRVSDNNtoASB <- harmonise_data(exposure_dat = exposure_HRV_SDNN_dat5e8_clumped,
                                 outcome_dat = outcome_ASB_dat, 
                                 action = 2)

table(DataMR_HRVSDNNtoASB$mr_keep) #TRUE 5

# This columns tells you which SNPs will be kept for the main analysis: All rows that are set to TRUE will be included in the MR analysis. 

attr(DataMR_HRVSDNNtoASB, "log") # Detailed summary of what was done and reasons for excluding SNPs 


# second threshold
DataMR_HRVSDNNtoASB_lib <- harmonise_data(exposure_dat = exposure_HRV_SDNN_dat5e5_clumped,
                                 outcome_dat = outcome_ASB_dat, 
                                 action = 2)

table(DataMR_HRVSDNNtoASB_lib$mr_keep) # FALSE 2 TRUE 52


# F statistics

DataMR_HRVSDNNtoASB$fstatistic <- ((DataMR_HRVSDNNtoASB$beta.exposure^2) / (DataMR_HRVSDNNtoASB$se.exposure^2))

nrow(DataMR_HRVSDNNtoASB[DataMR_HRVSDNNtoASB$fstatistic>10, ]) # 5 (so all!)
mean(DataMR_HRVSDNNtoASB$fstatistic) # 57.21474
summary(DataMR_HRVSDNNtoASB$fstatistic)
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #  33.19   50.61   51.77   57.21   53.87   96.64 


# liberal threshold

DataMR_HRVSDNNtoASB_lib <- DataMR_HRVSDNNtoASB_lib[DataMR_HRVSDNNtoASB_lib$mr_keep == TRUE, ]
DataMR_HRVSDNNtoASB_lib$fstatistic <- ((DataMR_HRVSDNNtoASB_lib$beta.exposure^2) / (DataMR_HRVSDNNtoASB_lib$se.exposure^2))

nrow(DataMR_HRVSDNNtoASB_lib[DataMR_HRVSDNNtoASB_lib$fstatistic>10, ]) # 52 (so all!)
mean(DataMR_HRVSDNNtoASB_lib$fstatistic) # 23.23841
summary(DataMR_HRVSDNNtoASB_lib$fstatistic)
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #  16.78   17.69   19.00   23.24   22.39   96.64  
```


``` r

# Perform standard MR analyses
(MR_HRVSDNNtoASB <- mr(DataMR_HRVSDNNtoASB, 
                  method_list = c("mr_weighted_median","mr_egger_regression",
                                  "mr_ivw")))

knitr::kable(as.data.frame(MR_HRVSDNNtoASB), "markdown")


# more liberal threshold
(MR_HRVSDNNtoASB_lib <- mr(DataMR_HRVSDNNtoASB_lib, 
                  method_list = c("mr_weighted_median", "mr_egger_regression",
                                  "mr_ivw")))

```


If you want to export your results in an excel file, you can run the following in R:
``` r
write.csv(MR_HRVSDNNtoASB, file = "output/Additional/MR_HRV_SDNNtoASB_unfiltered.csv") # Saves the results 
write.csv(MR_HRVSDNNtoASB_lib, file = "output/Additional/MR_HRV_SDNNtoASBlib_unfiltered.csv") # Saves the results 

```


### C. RMSSD

``` r
# Read HRV as exposure data

## Import the summary statistics file into R
SumStatsHRV <- fread(here::here("Data", 'HeartRateVarRMSSD_MR_Nolte_2017.txt') , header = T, data.table = T) 

head(SumStatsHRV)

# The summary statistics containing the converted effect size can now be exported and saved as a new file in our "data" folder
write.table(SumStatsHRV, "Data/HeartRateVarRMSSD_MR_Nolte_2017.txt", col.names = T, row.names = F, quote = F, sep = "\t") 

# read data using dedicated function of the TwoSampleMR package (pendent to the read_outcome_data() function)
exposure_HRV_RMSSD_dat <- read_exposure_data(filename = "Data/HeartRateVarRMSSD_MR_Nolte_2017.txt", 
                                      sep = "\t",
                                      snp_col = "SNP",           
                                      beta_col = "BETA",         
                                      se_col = "SE",       
                                      effect_allele_col = "A1",  
                                      other_allele_col = "A2",   
                                      pval_col = "P",  
                                      eaf_col = "MAF",
                                      samplesize_col = "N") 

# Provide a name for the outcome  
exposure_HRV_RMSSD_dat$exposure <- "HRV_RMSSD"     


# How many SNPs?
length(exposure_HRV_RMSSD_dat$SNP) # N = 2532340 SNPs
```


``` r}
# Select SNPs for exposure

# Select genome-wide significant SNPs
exposure_HRV_RMSSD_dat5e8 <- subset(exposure_HRV_RMSSD_dat, pval.exposure < 5e-8)


# Select a more liberal threshold
exposure_HRV_RMSSD_dat5e5 <- subset(exposure_HRV_RMSSD_dat, pval.exposure < 5e-5)

# Lets have a look at the number of SNPs that we would include:
length(exposure_HRV_RMSSD_dat5e8$SNP) # 205
length(exposure_HRV_RMSSD_dat5e5$SNP) # 1401

```



``` r}
# Keep all SNPs from the exposure set which are present in the outcome set
exposure_HRV_RMSSD_dat5e8_out <- exposure_HRV_RMSSD_dat5e8[exposure_HRV_RMSSD_dat5e8$SNP %in% outcome_ASB_dat$SNP, ]

exposure_HRV_RMSSD_dat5e5_out <- exposure_HRV_RMSSD_dat5e5[exposure_HRV_RMSSD_dat5e5$SNP %in% outcome_ASB_dat$SNP, ]


# Lets have a look at the number of SNPs we excluded
length(exposure_HRV_RMSSD_dat5e8_out$SNP) - length(exposure_HRV_RMSSD_dat5e5_out$SNP) # N excluded 541

```



``` r
# Clump the data

# using default (recommended) parameters
exposure_HRV_RMSSD_dat5e8_clumped <- clump_data(exposure_HRV_RMSSD_dat5e8_out)

# Lets have a look at the number of SNPs we excluded
length(exposure_HRV_RMSSD_dat5e8_clumped$SNP) # N included = 5
length(exposure_HRV_RMSSD_dat5e8_clumped$SNP) - length(exposure_HRV_RMSSD_dat5e8_out$SNP) # N excluded - 61


# second threshold

exposure_HRV_RMSSD_dat5e5_clumped <- clump_data(exposure_HRV_RMSSD_dat5e5_out)

# Lets have a look at the number of SNPs we excluded
length(exposure_HRV_RMSSD_dat5e5_clumped$SNP) # N included = 56

```


``` r
# Harmonise summary statistics
DataMR_HRVRMSSDtoASB <- harmonise_data(exposure_dat = exposure_HRV_RMSSD_dat5e8_clumped,
                                 outcome_dat = outcome_ASB_dat, 
                                 action = 2)

table(DataMR_HRVRMSSDtoASB$mr_keep) #TRUE 5

# This columns tells you which SNPs will be kept for the main analysis: All rows that are set to TRUE will be included in the MR analysis. 

attr(DataMR_HRVRMSSDtoASB, "log") # Detailed summary of what was done and reasons for excluding SNPs 


# second threshold
DataMR_HRVRMSSDtoASB_lib <- harmonise_data(exposure_dat = exposure_HRV_RMSSD_dat5e5_clumped,
                                 outcome_dat = outcome_ASB_dat, 
                                 action = 2)

table(DataMR_HRVRMSSDtoASB_lib$mr_keep) # FALSE 4 TRUE 52


# F statistics

DataMR_HRVRMSSDtoASB$fstatistic <- ((DataMR_HRVRMSSDtoASB$beta.exposure^2) / (DataMR_HRVRMSSDtoASB$se.exposure^2))

nrow(DataMR_HRVRMSSDtoASB[DataMR_HRVRMSSDtoASB$fstatistic>10, ]) # 5 (so all!)
mean(DataMR_HRVRMSSDtoASB$fstatistic) # 69.65482
summary(DataMR_HRVRMSSDtoASB$fstatistic)
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #  32.66   34.39   54.65   69.65   95.86  130.73


# liberal threshold

DataMR_HRVRMSSDtoASB_lib <- DataMR_HRVRMSSDtoASB_lib[DataMR_HRVRMSSDtoASB_lib$mr_keep == TRUE, ]
DataMR_HRVRMSSDtoASB_lib$fstatistic <- ((DataMR_HRVRMSSDtoASB_lib$beta.exposure^2) / (DataMR_HRVRMSSDtoASB_lib$se.exposure^2))

nrow(DataMR_HRVRMSSDtoASB_lib[DataMR_HRVRMSSDtoASB_lib$fstatistic>10, ]) # 52 (so all!)
mean(DataMR_HRVRMSSDtoASB_lib$fstatistic) # 24.48153
summary(DataMR_HRVRMSSDtoASB_lib$fstatistic)
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #  16.71   17.48   19.04   24.48   21.52  130.73  
```


``` r

# Perform standard MR analyses
(MR_HRVRMSSDtoASB <- mr(DataMR_HRVRMSSDtoASB, 
                  method_list = c("mr_weighted_median","mr_egger_regression",
                                  "mr_ivw")))

knitr::kable(as.data.frame(MR_HRVRMSSDtoASB), "markdown")


# more liberal threshold
(MR_HRVRMSSDtoASB_lib <- mr(DataMR_HRVRMSSDtoASB_lib, 
                  method_list = c("mr_weighted_median", "mr_egger_regression",
                                  "mr_ivw")))

```


If you want to export your results in an excel file, you can run the following in R:
``` r
write.csv(MR_HRVRMSSDtoASB, file = "output/Additional/MR_HRV_RMSSDtoASB_unfiltered.csv") # Saves the results 
write.csv(MR_HRVRMSSDtoASB_lib, file = "output/Additional/MR_HRV_RMSSDtoASBlib_unfiltered.csv") # Saves the results 

```


## Multivariate MR analyses with RHR and HRV as exposures


``` r}
# Read in RHR exposure
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

# Select genome-wide significant SNPs
exposure_HR_dat5e8 <- subset(exposure_HR_dat, pval.exposure < 5e-8)

# Keep all SNPs from the exposure set which are present in the outcome set
exposure_HR_dat5e8_out <- exposure_HR_dat5e8[exposure_HR_dat5e8$SNP %in% outcome_ASB_dat$SNP, ]

# Use HRV RMSSD exposure from above
exposure_HRV_RMSSD_dat5e8_out

```



``` r}
## combine the two single data sets to get SNPs that are significant for at least one exposure
combined_snps <- unique(rbind(subset(exposure_HR_dat5e8_out, select = c(SNP)), subset(exposure_HRV_RMSSD_dat5e8_out, select = c(SNP))))

rownames(combined_snps) <- NULL

## filter two original data sets to only get significant SNPs (for at least one exposure)
exposure_HR_dat_sig <- exposure_HR_dat %>% subset(exposure_HR_dat$SNP %in% combined_snps$SNP)

exposure_HRV_dat_sig <- exposure_HRV_RMSSD_dat %>% subset(exposure_HRV_RMSSD_dat$SNP %in% combined_snps$SNP)


# clump again to exclude correlated SNPs among different exposure SNPs
exposure_HR_dat_clumped <- clump_data(exposure_HR_dat_sig, clump_r2 = 0.001, clump_p1 = 1, clump_p2 = 1, clump_kb = 10000)

exposure_HRV_dat_clumped <- clump_data(exposure_HRV_dat_sig, clump_r2 = 0.001, clump_p1 = 1, clump_p2 = 1, clump_kb = 10000)

# Remove sample size for exposure from HRV
exposure_HRV_dat_clumped <- exposure_HRV_dat_clumped %>% 
    select(-samplesize.exposure)

## combine two data sets with final SNPs
exposure_HRmeasures_clumped <- rbind(exposure_HR_dat_clumped, exposure_HRV_dat_clumped)


## filter out one SNP with wrong direction of effect  (RHR on ASB) (from previous Steiger filtering)

wrongDirect <- c("rs1015149", "rs10895277", "rs11229113", "rs11788032", "rs12174018",
                 "rs13129601", "rs16835702", "rs2162889", "rs2246363", "rs2586886",
                 "rs6457796", "rs7193221", "rs73135307", "rs7795919") 



exposure_HRmeasures_clumped <- exposure_HRmeasures_clumped[!(exposure_HRmeasures_clumped$SNP %in% wrongDirect), ]

## harmonize with outcome data
mvDat_HR_clumped <- mv_harmonise_data(exposure_dat = exposure_HRmeasures_clumped, outcome_dat = outcome_ASB_dat, harmonise_strictness = 2)


##############


# perform MR
(MultiMr_clumped <- mv_multiple(mvdat = mvDat_HR_clumped, intercept = FALSE, instrument_specific = FALSE, pval_threshold = 1))

write.csv(MultiMr_clumped, file = "output/Additional/MultiMr_clumped.csv") # Saves the results in dataframe

```




## Univariable MR analyses between RHR and HRV (positive control)

``` r
## read outcome data

# HRV
outcome_HRV_dat <- read_outcome_data(snps = NULL,                            # Keep 'snps = NULL' so that all SNPs are uploaded
                                     filename = "Data/HeartRateVarRMSSD_MR_Nolte_2017.txt",
                                     sep = "\t",                             # file tab ("\t") or  space delimited ("")
                                     snp_col = "SNP",                        # Insert the name of column that includes the SNP IDs
                                     beta_col = "BETA",                      # Insert the name of column with effect sizes
                                     se_col = "SE",                          # Insert the name of column with standard errors
                                     effect_allele_col = "A1",               # Insert the name of column with the effect allele
                                     other_allele_col = "A2",                # Insert the name of column with the non-effect allele
                                     pval_col = "P",                         # Insert the name of column with p-value
                                     eaf_col = "MAF",                        # effect allele frequency
                                     samplesize_col = "N")                   # sample size available for each SNP

# Note: A few columns are not present in the sum stat file so cannot be read directly but can be added

# Provide a name for the outcome
outcome_HRV_dat$outcome <- "HRV" 

head(outcome_HRV_dat)

# How many SNPs are included?
length(outcome_HRV_dat$SNP) # Number of SNPs included 9590567

# Check the number of rows with command line
system("wc -l Data/HeartRateVarRMSSD_MR_Nolte_2017.txt") # one more than the R object as previously as counts column names as one row 
system("head Data/HeartRateVarRMSSD_MR_Nolte_2017.txt")
system("tail Data/HeartRateVarRMSSD_MR_Nolte_2017.txt")

```




``` r}
# Keep all SNPs from the exposure set which are present in the outcome set
exposure_HR_dat5e8_out <- exposure_HR_dat5e8[exposure_HR_dat5e8$SNP %in% outcome_HRV_dat$SNP, ]

# Lets have a look at the number of SNPs we excluded
length(exposure_HR_dat5e8_out$SNP) - length(exposure_HR_dat5e8$SNP) # N excluded -31634

```



``` r
# ========= Clump the data =========

# using default (recommended) parameters
exposure_HR_dat5e8_clumped <- clump_data(exposure_HR_dat5e8_out)

# Lets have a look at the number of SNPs we excluded
length(exposure_HR_dat5e8_clumped$SNP) # N included = 241
length(exposure_HR_dat5e8_clumped$SNP) - length(exposure_HR_dat5e8_out$SNP) # N excluded - 14859


# F statistic
# ((n-k-1)/k) * (R2/ 1-R2)
# ((458835-297-1)/297) * (0.01 / (1 - 0.01))
# ((458835-297-1)/297) * (0.02 / (1 - 0.02)) 
# ((458835-297-1)/297) * (0.05 / (1 - 0.05)) 
```



``` r
####### Harmonize
DataMR_HRtoHRV <- harmonise_data(exposure_dat = exposure_HR_dat5e8_clumped,
                                 outcome_dat = outcome_HRV_dat, 
                                 action = 2)

table(DataMR_HRtoHRV$mr_keep)  # FALSE 12 TRUE 261
# This columns tells you which SNPs will be kept for the main analysis: All rows that are set to TRUE will be included in the MR analysis. 

attr(DataMR_HRtoHRV, "log") # Detailed summary of what was done and reasons for excluding SNPs 

# F statistics

DataMR_HRtoHRV$fstatistic <- ((DataMR_HRtoHRV$beta.exposure^2) / (DataMR_HRtoHRV$se.exposure^2))

nrow(DataMR_HRtoHRV[DataMR_HRtoHRV$fstatistic>10, ]) # 273 
mean(DataMR_HRtoHRV$fstatistic) # 84.81837
summary(DataMR_HRtoHRV$fstatistic)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  26.84   36.82   49.03   84.82   77.21 1185.15  
```

 
``` r

# Run as subset of estimators
(MR_HRtoHRV <- mr(DataMR_HRtoHRV, 
                  method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression")))

knitr::kable(as.data.frame(MR_HRtoHRV), "markdown")


# mr_raps

mr_frame1 <- as.data.frame(DataMR_HRtoHRV[DataMR_HRtoHRV$mr_keep == TRUE, # keep only SNPS after harmoization
                                          c('beta.exposure', 'se.exposure', 'beta.outcome', 'se.outcome')])
names(mr_frame1) <- c("RHR","RHR_se","HRV","HRV_se")

# Perform MR-RAPS 
(mr_raps_result <- mr.raps.overdispersed.robust(b_exp = mr_frame1$RHR, b_out = mr_frame1$HRV, # define effects of SNPs on exposure and outcome
                                                se_exp = mr_frame1$RHR_se, se_out = mr_frame1$HRV_se, # define standard errors of exposure and outcome
                                                loss.function = "huber", # choose loss function from c("huber", "tukey")
                                                k = 1.345, # Threshold parameter in the Huber and Tukey loss functions. (default)
                                                initialization = c("l2"), # Method to initialize the robust estimator, c("l2", "mode")
                                                suppress.warning = FALSE, 
                                                niter = 20, # Maximum number of interations to solve the estimating equations. (default)
                                                tol = .Machine$double.eps^0.5)) # Numerical precision. (default)

knitr::kable(as.data.frame(mr_raps_result), "markdown")



```


If you want to export your results in an excel file, you can run the following in R:

``` r
write.csv(MR_HRtoHRV, file = "output/Controls/MR_restingHRtoHRV.csv") # Save the results 

```


