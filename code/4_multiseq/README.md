# 705 Cox Multi-Sequence Sieve Analysis

This document implements a multi-sequence sieve analysis using methods from  
Peng et al. (in progress) to analyze binary within-host sequence marks.

## Data Notes
- 2,188 participants total  
- 119 with `hiv1event == 1`  
- 23 controls with sequencing data (`hiv1event == 0`)  
- 3 cases without sequencing data (removed)  
- Final analysis dataset: 2,185 participants  

Covariates included: **BMI**, **risk score**, **age**, with stratification by **South Africa indicator**.

## Workflow

### 1. Data Processing
- Reads the multi-sequence, mindist, and subject-level sieve datasets  
- Restricts to the first visit per subject  

### 2. Pre-Screening of Binary Marks
For each binary mark:
- Compute raw within-host proportion \( Q = K / \text{sequence\_depth} \)  
- Create thresholded indicators \( J_{q_0} \) at **q₀ = 0.01** and **0.99**  
- Compare thresholded marks to mindist marks  
- Screen in a mark if:
  1. ≥12 disagreements with mindist marks among cases  
  2. ≥4 cases with \( J_{q_0} = 1 \)  
  3. ≥4 cases with \( J_{q_0} = 0 \)

Pre-screening outputs CSV files for:
- Insert-matched marks  
- Reference-independent AA scanning marks  
- Sequon marks  

### 3. Cox Sieve Analysis
For each screened-in mark, the script runs `run_full_analysis()`:
- Covariates: treatment (`Z`), age, BMI, risk score  
- Event/time: `allstatus`, `observed_time`  
- Stratification by South Africa  
- Bootstrap SEs (default: 300 samples)

Outputs for each mark:
- VE estimates for \( J_{q_0} = 0 \) and 1  
- Bootstrap SEs  
- Confidence intervals  
- P-values for:  
  - VE₀ = 0  
  - VE₁ = 0  
  - VE₀ = VE₁ = 0  
  - VE₀ = VE₁  

CSV files produced:
- `res_match_0.01.csv`, `res_match_0.99.csv`  
- `res_aa_0.01.csv`, `res_aa_0.99.csv`

### 4. Sequon Marks
Sequon pre-screening is run but typically yields **no screened-in features** at either threshold.

### 5. Permutation Test 
A slow but complete permutation procedure is provided to compute familywise error–adjusted p-values.  
It permutes case sequence depths and counts, re-runs the sieve model for all marks, and applies `p.adj.perm2()`.

### 6. Significant mark: hxb2.832.is.L
The document includes:
- Code for plotting raw within-host prevalence of L832  
- Sensitivity analysis excluding gaps  
- Export of graph PDFs

For any questions, please contact James Peng (jpspeng@uw.edu). 