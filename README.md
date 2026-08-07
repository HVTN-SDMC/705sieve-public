## Genotypic sieve analysis in the HVTN 705/HPX2008 Imbokodo trial
*R* code implementing data analyses that generated figures and tables in the manuscript

Juraska, Li et al., *Quantifying how HIV-1 Envelope Sequence Features Impact Vaccine Efficacy in the HVTN 705/HPX2008 Randomized Trial in Southern African Women*.

### 1. System requirements and installation guide 

Unless otherwise specified, the software was tested on macOS Sequoia 15.6.1 running *R* version 4.4.3.
R packages required for analyses are listed in the renv.lock file. To restore all the *R* packages specified in the renv.lock file, follow the instructions below:

 *  Install required version of *R*.  
 * Download the repository to your local machine and start a project in the Rstudio in the folder.
 * Install *R* package "renv".  
 * Run renv::restore(). It creates the folder *renv/library/* and installs packages or links recorded in renv.lock file into the project library.
  
  
### 2 Main programs and *R* files required for each analysis

* `code/0_preScreen/screenMain.R`: performs the minimum variability filter  
* `code/1_westfallYoung/sieveWestfallYoungPermPvals.R`: performs sieve analyses and p-value adjustment
* `code/2_sieveBinary/sieveBinaryMain.R`: summarizes sieve analyses for binary features of Env amino acid sequences
* `code/3_sievePH/VEbyHammingDist.R` and `code/3_sievePH/VEbyOtherQuantMarks.R`: plot sieve analysis results for continuous-valued features of Env amino acid sequences
* `code/4_multiseq/multiseq_heatmaps.R`: generates heatmaps based on multiple sequences per participant
* `code/4_multiseq/cox_multiseq_hvtn705_sieve_analysis.Rmd`: performs multi-sequence sieve analysis
* `code/5_viralLoad/viralLoad.R`: generates violin and boxplots for viral loads and sequencing depths and computes p-values for comparing vaccine to placebo
* `code/6_TMscore/TMscoreWestfallYoungPermPvals.R`: TM score analysis comparing the TM scores of vaccine endpoints to placebo endpoints using linear regressions and p-value adjustment
* `code/6_TMscore/TMscore.R`: plots TM score analysis 
* `code/7_diproperm/1_data_prep_dpp.R` and `code/7_diproperm/2_analyze_dpp.R`: performs the Diproperm high-dimensional feature analysis

### 3. How to run the *R* files

  Each script file should run in less than 10 minutes except for `code/1_westfallYoung/sieveWestfallYoungPermPvals.R` and 
  `code/6_TMscore/TMscoreWestfallYoungPermPvalse.R` which may require several hours to complete.  Outputs are saved
  as either figures (PDF files in the `figures` directory) or tables (CSV files in the `tables` directory). 
  Outputs from `code/4_multiseq/cox_multiseq_hvtn705_sieve_analysis.Rmd` are the same folder as the code.
  
  From the command line, starting in the `code` directory, run the following commands:

    R CMD BATCH 0_preScreen/screenMain.R &
    R CMD BATCH 1_westfallYoung/sieveWestfallYoungPermPvals.R &
    R CMD BATCH 2_sieveBinary/sieveBinaryMain.R &
    R CMD BATCH 3_sievePH/VEbyHammingDist.R &
    R CMD BATCH 3_sievePH/VEbyOtherQuantMarks.R &
    R CMD BATCH 4_multiseq/multiseq_heatmaps.R &
    Rscript -e "rmarkdown::render('4_multiseq/cox_multiseq_hvtn705_sieve_analysis.Rmd')" &
    R CMD BATCH 5_viralLoad/viralLoad.R &
    R CMD BATCH 6_TMscore/TMscoreWestfallYoungPermPvalse.R &
    R CMD BATCH 6_TMscore/TMscore.R
    R CMD BATCH 7_diproperm/1_data_prep_dpp.R
    R CMD BATCH 7_diproperm/2_analyze_dpp.R

### 4. Mapping from Code to main figures and tables
#### Table 3
 - `code/0_preScreen/screenMain.R` generates the table `tables/0_preScreen/posIsAAtier1posScreenedIn.csv` which contains the Tier 1 screened-in residue presence/absence features for sieve analysis.  

#### Figure 1
 - `code/5_viralLoad/viralLoad.R` generates the plots in all the panels  
 - `figures/5_viralLoad/viralLoad_trend_lmm.pdf` which contains the longitudinal trajectories of HIV-1 viral load (copies/ml) in primary endpoints in the first RNA-positive sample and 2–6 weeks later while pre-ART initiation
 - `figures/5_viralLoad/vl_ave_vlbyTreatment.pdf` which contains the boxplots of HIV-1 viral load (copies/ml) in primary endpoints
 - `figures/5_viralLoad/NTSeqdepthbyTreatment.pdf` which contains the boxplots of Env sequencing depth for nucleotide sequences in primary endpoints
 - `figures/5_viralLoad/AASeqdepthbyTreatment.pdf` which contains the boxplots of Env sequencing depth for amino acid sequences in primary endpoints

#### Figure 2
 - `code/1_westfallYoung/sieveWestfallYoungPermPvals.R` generates the adjusted p-values in the panel A
   - `tables/westfallYoung/WestfallYoungAdjPvalues_tier1Type1to4.csv` contains the adjusted p-values 
 - `code/2_sieveBinary/sieveBinaryMain.R` generates the vaccine efficacy (VE) by amino acid residue at each Env position in the panel A and plots the forest plot for Env 364 in panel C
   - `tables/sieveBinary/tier1posIsAAVE_subTypeE.csv` contains the VE by screened-in amino acid residue at each Env position in Tier 1 hypothesis-driven Env amino acid positions
 - `code/2_sieveBinary/volcanoPlot.R` generates the volcano plot in the panel A
   - `figures/sieveBinary/tier1_posIs_Volcano.pdf` contains the volcano plot 

#### Figure 3 
 - `code/1_westfallYoung/sieveWestfallYoungPermPvals.R` generates the adjusted p-values in both panels
   - `tables/westfallYoung/WestfallYoungAdjPvalues_tier1Type5to7.csv` contains the adjusted p-values 
 - `code/3_sievePH/VEbyHammingDist.R` generates the plots in both panels
   - `figures/sievePH/705_sievePH_VE_hdist_zspace_c97za_c_ab_PP.pdf` contains the VE against viruses with physicochemical (PC)-weighted Hamming distances from the C97ZA Env vaccine insert sequence in clade C bnAb resistance-associated signature positions 
   - `figures/sievePH/705_sievePH_VE_hdist_zspace_c97za_hvtn505_cd4bs_kmer_PP.pdf` contains the VE against viruses with physicochemical (PC)-weighted Hamming distances from the C97ZA Env vaccine insert sequence in CD4 binding site-overlapping HVTN 505 signature k-mers
   
#### Figure 4
 - `code/4_multiseq/cox_multiseq_hvtn705_sieve_analysis.Rmd` generates the figures in both panels
   - `code/4_multiseq/res_aa_0.99.csv`contains the VE against viral populations with $\geq 99%$ vs. $<99%$ prevalence of leucine at Env position 832 and VE against viral populations whose single representative sequence carries L832 vs. notL832
   - `code/4_multiseq/L832_prevalence_plot_gaps_included.pdf` contains the observed L832 prevalence among sampled sequences from primary endpoints by treatment group
   - ? contains the adjusted P-values 
   