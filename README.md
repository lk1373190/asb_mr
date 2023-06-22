
# Resting heart rate and antisocial behaviour: A Mendelian randomisation study.


## Description
The current repository contains scripts and data for reproducing the findings from the article ["Resting heart rate and antisocial behaviour: A Mendelian randomisation study."](INSERT LINK HERE). 

## Scripts
Scripts are provided in .md formats.
* Main analyses - "main_hr_asb_mr_script.md"
* Additional analyses - "additional_hr_asb_mr_script.md

### Package info
The following packages (and versions) are required to run the code:
* arrangements (1.1.9)
* data.table (1.14.8)
* devtools (2.4.5)
* dpylr (1.1.2)
* gmodels (2.18.1.1)
* here (1.0.1)
* MendelianRandomization (0.7.0)
* mr.raps (0.2)
* MRPRESSO (1.0)
* TwoSampleMR (0.5.7)

### Sources
A lot of the code used in this project has been inspired by the "TwoSampleMR" package and tutorial by
[@mrcieu](https://mrcieu.github.io/TwoSampleMR/articles/introduction.html)

## Data
No data files are included in this repository. The data analysed in the MR analysis were sourced from the following Genome Wide Association Studies.
### 1. Antisocial behaviour
* [Uncovering the genetic architecture of broad antisocial behavior through a genome-wide association study meta-analysis](https://doi.org/10.1038/s41380-022-01793-3)
* Data available on request directly from the study authors: [Dr Jorim J. Tielbeek](j.tielbeek@amsterdamumc.nl).

### 2. Resting heart rate
* [Genetic overlap of chronic obstructive pulmonary disease and cardiovascular disease-related traits: a large-scale genomewide cross-trait analysis](https://doi.org/10.1186/s12931-019-1036-8)
* Data (controlling for smoking) available on request from directly from the study authors [Dr Zhaozhong Zhu](zhz586@mail.harvard.edu)

### 3. Heart rate variability
* [Genetic loci associated with heart rate variability and their effects on cardiac disease risk)[https://doi.org/10.1038/ncomms15805]
* Data available on from the [GWAS catalog](https://www.ebi.ac.uk/gwas/publications/28613276)
