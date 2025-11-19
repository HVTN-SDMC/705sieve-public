## Genotypic Sieve Analysis in the HVTN 705/HPX2008 Imbokodo Trial
*R* code implementing data analyses that generated figures and tables in the manuscript Juraska, Li et al., Quantifying how a Mosaic HIV-1 Vaccine Regimen’s Efficacy Depends on Envelope Sequence Features in the HVTN 705/HPX2008 Imbokodo Trial.

### 1. System Requirements

  Unless otherwise specified, software was tested on macOS Monterey running R version 4.4.1 (2024-06-14).
  
  Main programs as well as specific R packages and R files required by each program are listed below.

* `code/0_preScreen/screenMain.R`: performs the minimum variablity filter  
    + tidyverse 2.0.0 
    + plyr 1.8.9 
    + screenUtils.R 
    + common.R 
    
* `code/1_westfallYoung/sieveWestfallYoungPermPvals.R`: performs sieve analyses and p-value adjustment
    + sievePH 1.0.3
    + lunnMcneil.R
    + p.adj.perm2.R
    + common.R
    
* `code/2_sieveBinary/sieveBinaryMain.R`: summarizes sieve analyses for binary features of Env amino acid sequences
    + tidyverse 2.0.0
    + plyr 1.8.9
    + grid 4.4.1
    + gridExtra 2.3
    + gtable 0.3.5
    + lunnMcneil.R
    + sieveBinaryUtils.R
    + forest.R
    + common.R
    
* `code/3_sievePH/VEbyHammingDist.R` and `code/3_sievePH/VEbyOtherQuantMarks.R`: plot sieve analysis results for continuous-valued features of Env amino acid sequences
    + tidyverse 2.0.0
    + sievePH 1.0.3
    + ggplot.summary.sievePH.R
    + common.R

* `code/4_multiseq/multiseq_heatmaps.R`: generates heatmaps based on multiple sequences per participant
    + tidyverse 2.0.0
    + ggpubr 0.6.0
    + scales 1.3.0
    + DT 0.33.3
* `code/4_multiseq/cox_multiseq_hvtn705_sieve_analysis.Rmd`: performs multi-sequence sieve analysis
    + deconv_fast.R
    + multiseq_cox_code.R
    + p.adj.perm2.R
    + dplyr 1.1.4
    + DT 0.33.3
    + purrr 1.0.4
    + ggplot2 4.0.0
    + ggpubr 0.6.0
    + survival 3.8-3
    + deconvolveR 1.2-1
    + future 1.58.0
    + future.apply 1.20.0
    + ggbeeswarm 0.7.2
    
* `code/5_viralLoad/viralLoad.R`: generates violin and boxplots for viral loads and sequencing depths and computes p-values for comparing vaccine to placebo
    + tidyverse 2.0.0
    + scales 1.3.0

* `code/6_TMscore/TMscoreWestfallYoungPermPvalse.R`: TM score analysis comparing the TM scores of vaccine endpoints to placebo endpoints using linear regressions and p-value adjustment
    + tidyverse 2.0.0

* `code/6_TMscore/TMscore.R`: plot TM score analysis 
    + tidyverse 2.0.0

* `code/7_diproperm/1_data_prep_dpp.R` and `code/7_diproperm/2_analyze_dpp.R`: performs Diproperm high-dimensional feature analysis
    + https://github.com/youyifong/diproperm

### 2. Installation Guide
  
* Install required version of *R*.  
* Install required *R* packages.  
* Clone this repository.
  
### 3. User Instructions

  Each code file should run in less than 10 minutes except `code/1_westfallYoung/sieveWestfallYoungPermPvals.R` and 
  `code/6_TMscore/TMscoreWestfallYoungPermPvalse.R` which may take approximately a few hours to run.  Output is saved
  as either a figure (a `.pdf` file in the `figures` directory) or a table (a `.csv` file in the `tables` directory).
  
  From the command line, starting in the `code` directory, run the following commands:

    R CMD BATCH 0_preScreen/screenMain.R &
    R CMD BATCH 1_westfallYoung/sieveWestfallYoungPermPvals.R &
    R CMD BATCH 2_sieveBinary/sieveBinaryMain.R &
    R CMD BATCH 3_sievePH/VEbyHammingDist.R &
    R CMD BATCH 3_sievePH/VEbyOtherQuantMarks.R &
    R CMD BATCH 4_multiseq/multiseq_heatmaps.R &
    R CMD BATCH 4_multiseq/cox_multiseq_hvtn705_sieve_analysis.Rmd &
    R CMD BATCH 5_viralLoad/viralLoad.R &
    R CMD BATCH 6_TMscore/TMscoreWestfallYoungPermPvalse.R &
    R CMD BATCH 6_TMscore/TMscore.R
    R CMD BATCH 7_diproperm/1_data_prep_dpp.R
    R CMD BATCH 7_diproperm/2_analyze_dpp.R
    
    
