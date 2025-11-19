

# ---------------------------------------------------------------------------- #
# STEP 0:  set up our working environment
# ---------------------------------------------------------------------------- #


# refresh workspace
rm (list=ls ())
set.seed (378)

# load required packages
#library ("diproperm")


# ---------------------------------------------------------------------------- #
# STEP 1:  load and filter our indata
# ---------------------------------------------------------------------------- #


# where does our data live?
path.home <- "/path/to/home/directory"

# load ptids of primary endpoints to include in analysis
data.ptid <- read.csv (file.path (path.home, "ptid.csv"), header=T)[, 1]

# load marks file and filter by primary endpoints
data.marks <- read.csv (file.path (path.home, "marks/HVTN_705_sieve_data_final_v4.csv"), header=T)
data.marks <- data.marks[data.marks$subjid %in% data.ptid, ]

# load structural/AF3 marks:  windowed TM scores
data.tm <- read.csv (file.path (path.home, "tm_score/tm_score_window_marks_lb6_v705_A_v1.csv"), header=T)

# load structural/AF3 marks:  region-based TM scores
data.tm.region <- read.csv (file.path (path.home, "tm_score/tm_score_region_marks_v705_A_v1.csv"), header=T)

# merge our TM scores with our marks data into an analysis-ready dataset
data.analysis <- merge (data.marks, data.tm, by.x="subjid", by.y="ptid")
data.analysis <- merge (data.analysis, data.tm.region, by.x="subjid", by.y="ptid")


# ---------------------------------------------------------------------------- #
# STEP 2:  subset our data for analysis
# ---------------------------------------------------------------------------- #


#  SET 0. Treatment assignments (the class or response variable)
data.treatment <- as.numeric (data.analysis$armdesc == "Vaccine")

#  SET 1. Control variables only
#    -- Age, BMI, baseline risk score and enrollment in South Africa
data.set.1 <- data.analysis[, c ("age", "bmi", "riskscoresl", "ind_sa")]

#  SET 2. Tier 1:  Env gp120 C97ZA insert residue match/mismatch status
#    (a) V1V2-hotspot Janssen NHP challenge: sites 157-184
data.set.2a <- data.analysis[, 2331:2361]

#    (b) IgG3 V1V2 breadth-correlates sites: sites 130, 155, 160, 161, 164, 165, 166, 170, 172, 173, 178, 179, 181, 186, 200
data.set.2b <- data.analysis[, c (2326, 2329, 2334, 2335, 2338:2340, 2344, 2346, 2349, 2355, 2356, 2358, 2379, 2388)]

#    (c) RV144 V2 C1 region sites within ADCC epitope: sites 51-61, 69-78
data.set.2c <- data.analysis[, c (2116:2127, 2135:2144)]

#    (d) RV144 primary sieve sites (Rolland, Edlefsen et al. 2012): sites 169 and 181
data.set.2d <- data.analysis[, c ("hxb2.169.1mer.c97za.match.tier1", "hxb2.181.1mer.c97za.match.tier1")]

#    (e) Signature sites for C-clade bnAbs (i.e., VRC01, 4E10, 3BNC117, and CAP256.VRC26):  sites 230, 279-282, 350, 364, 429, 432, 456, 471
data.set.2e <- data.analysis[, c (2420, 2501:2504, 2576, 2598, 2664, 2667, 2694, 2721)]

#  SET 3. Tier 1:  Env gp70-1428 V1V2 antigen sequence residue match/mismatch status
#    -- gG3 V1V2 breadth correlates sites: sites 130, 155, 160, 161, 164, 165, 166, 170, 172, 173, 178, 179, 181, 186, 200
data.set.3 <- data.analysis[, 3052:3066]

#  SET 4. Tier 1:  Env gp120 C97ZA vaccine insert-independent residue presence/absence indicators
#    (a) V1V2-hotspot Janssen NHP challenge: sites 157-184
data.set.4a <- data.analysis[, 3721:3890]

#    (b) IgG3 V1V2 breadth-correlates sites: sites 130, 155, 160, 161, 164, 165, 166, 170, 172, 173, 178, 179, 181, 186, 200
data.set.4b <- data.analysis[, c (3678:3689, 3704:3712, 3726:3737, 3747:3767, 3786:3792, 3803:3812, 3817:3827, 3841:3859, 3861:3865, 4024:4038, 4064:4069)]

#    (c) RV144 V2 C1 region sites within ADCC epitope: sites 51-61, 69-78
data.set.4c <- data.analysis[, c (3331:3353, 3379:3395)]

#    (d) RV144 primary sieve sites (Rolland, Edlefsen et al. 2012): sites 169 and 181
data.set.4d <- data.analysis[, c (3776:3785, 3861:3865)]

#    (e) Signature sites for C-clade bnAbs (i.e., VRC01, 4E10, 3BNC117, and CAP256.VRC26):  sites 230, 279-282, 350, 364, 429, 432, 456, 471
data.set.4e <- data.analysis[, c (4136:4139, 4335:4362, 4739:4752, 4895:4899, 5057:5063, 5068:5070, 5158:5164, 5278:5285)]

#  SET 5. Tier 1:  Presence/absence of sequons
#    (a) Glycan holes in C97ZA insert (Korber CHAVD presentation, Oct. 2021): sites 234, 392
data.set.5a <- data.analysis[, c (6956, 6960)]

#    (b) V2 bnAb contact glycans: sites 130, 156, 160, 234, 332, 362, 363, 392
data.set.5b <- data.analysis[, 6953:6960]

#  SET 6. Tier 1:  PC-weighted Hamming distance from vaccine insert
#    (a) V1V2-hotspot Janssen NHP challenge: sites 157-184
data.set.6a <- data.frame (hdist.zspace.c97za.v2=data.analysis$hdist.zspace.c97za.v2)

#    (b) IgG3 V1V2 breadth correlates sites: sites 130, 155, 160, 161, 164, 165, 166, 170, 172, 173, 178, 179, 181, 186, 200
data.set.6b <- data.frame (hdist.zspace.1428v1v2=data.analysis$hdist.zspace.1428v1v2)

#    (c) RV144 V2 C1 region sites within ADCC epitope: sites 51-61 and 69-78
data.set.6c <- data.frame (hdist.zspace.c97za.adcc=data.analysis$hdist.zspace.c97za.adcc)

#    (d) Signature sites for C-clade bnAbs (i.e., VRC01, 4E10, 3BNC117, and CAP256.VRC26):  sites 230, 234, 279-282, 350, 364, 429, 432, 456, 471
data.set.6d <- data.frame (hdist.zspace.c97za.c_ab=data.analysis$hdist.zspace.c97za.c_ab)

#    (e) HVTN 505 sieve analysis (DeCamp et al.): 93 CD4bs mAb contact sites
data.set.6e <- data.frame (hdist.zspace.c97za.hvtn505.cd4bs.antibody=data.analysis$hdist.zspace.c97za.hvtn505.cd4bs.antibody)

#    (f) HVTN 505 sieve analysis (DeCamp et al.): 54 sites corresponding to 9-mers in 4 linear regions overlapping the CD4bs contact residues
data.set.6f <- data.frame (hdist.zspace.c97za.hvtn505.cd4bs.kmer=data.analysis$hdist.zspace.c97za.hvtn505.cd4bs.kmer)

#  SET 7. Tier 1:  PC-weighted Hamming dist from gp70-1428 V1V2 antigen
#    -- IgG3 V1V2 breadth correlates sites: sites 130, 155, 160, 161, 164, 165, 166, 170, 172, 173, 178, 179, 181, 186, 200
data.set.7 <- data.frame (hdist.zspace.1428v1v2=data.analysis$hdist.zspace.1428v1v2)

#  SET 8. Tier 1:  Other Env AA features
#    (a) Length of V1V2
data.set.8a <- data.frame (length.v1v2=data.analysis$length.v1v2)

#    (b) Number of sequons in V1V2
data.set.8b <- data.frame (num.sequons.v1v2=data.analysis$num.sequons.v1v2)

#    (c) Electrochemical charge of V2
data.set.8c <- data.frame (charge.v2=data.analysis$charge.v2)

#    (d) Cysteine count in gp160
data.set.8d <- data.frame (cysteine.count=data.analysis$cysteine.count)

#    (e) Indicators of single vs multiple viral lineages
data.set.8e <- data.frame (transmitted.founder.status=data.analysis$transmitted.founder.status)

#  SET 9. Tier 2:  Vaccine residue match/mismatch indicators in vaccine-insert region of Env, excluding Tier 1 hypothesis-driven position set
#    -- All sufficiently-variable positions in signal peptide + gp120 + N-terminal portion of gp41 included in vaccine inserts, excluding Tier 1 hypothesis-driven positions to avoid duplication
data.set.9 <- data.analysis[, 2050:3051]
data.set.9 <- data.set.9[, colSums (data.set.9) >= 4 & colSums (data.set.9) <= (nrow (data.set.9) - 3)]
data.set.9 <- data.set.9[, names (data.set.9)[!(names (data.set.9) %in% unique (c (names (data.set.2a), names (data.set.2b), names (data.set.2c), names (data.set.2d), names (data.set.2e))))]]

#  SET 10. Tier 2:  Vaccine insert-independent residue presence/absence indicators in all of Env, excluding Tier 1 hypothesis-driven position set
#    -- All sufficiently-variable positions in Env, excluding Tier 1 hypothesis-driven positions to avoid duplication (vaccine insert residues are ignored in the vaccine-insert region)
data.set.10 <- data.analysis[, 3067:6918]
data.set.10 <- data.set.10[, colSums (data.set.10) >= 4 & colSums (data.set.10) <= (nrow (data.set.10) - 3)]
data.set.10 <- data.set.10[, names (data.set.10)[!(names (data.set.10) %in% unique (c (names (data.set.4a), names (data.set.4b), names (data.set.4c), names (data.set.4d), names (data.set.4e))))]]

#  SET 11. Tier 2: PC-weighted Hamming distance to the vaccine insert
#    (a) All of gp120
data.set.11a <- data.frame (hdist.zspace.c97za.gp120=data.analysis$hdist.zspace.c97za.gp120)

#    (b) N-terminal portion of gp41 included in vaccine inserts
data.set.11b <- data.frame (hdist.zspace.c97za.gp41=data.analysis$hdist.zspace.c97za.gp41)

#    (c) All of V5
data.set.11c <- data.frame (hdist.zspace.c97za.v5=data.analysis$hdist.zspace.c97za.v5)

#  SET 12. Tier 2: Other AA features
#    (a) Length of V5
data.set.12a <- data.frame (length.v5=data.analysis$length.v5)

#    (b) Number of sequons in V5 (dichotomized as < 2 vs. >= 2)
data.set.12b <- data.frame (num.sequons.v5=as.numeric (data.analysis$num.sequons.v5 >= 2))

#  SET 13. Structural Features: TM score, comparing sequence with C97ZA insert
#    (a) All 15mer window TM scores
data.set.13a <- data.analysis[, 6966:7080]

#    (b) C1 region
data.set.13b <- data.frame (tmscore.region.C1=data.analysis$tmscore.region.C1)

#    (c) V1/V2 loop
data.set.13c <- data.frame (tmscore.region.V1.V2=data.analysis$tmscore.region.V1.V2)

#    (d) Loop D
data.set.13d <- data.frame (tmscore.region.loop.D=data.analysis$tmscore.region.loop.D)

#    (e) C2 region
data.set.13e <- data.frame (tmscore.region.C2=data.analysis$tmscore.region.C2)

#    (f) C3 region
data.set.13f <- data.frame (tmscore.region.C3=data.analysis$tmscore.region.C3)

#    (g) C4 region
data.set.13g <- data.frame (tmscore.region.C4=data.analysis$tmscore.region.C4)

#    (h) C5 region
data.set.13h <- data.frame (tmscore.region.C5=data.analysis$tmscore.region.C5)

#    (i) Items (b) through (g) combined
data.set.13i <- data.analysis[, c (7082:7087, 7093)]

#    (j) Beta2 strand (from the N-terminal V1/V2 stem): positions 121-123
data.set.13j <- data.frame (tmscore.region.beta2=data.analysis$tmscore.region.beta2)

#    (k) Beta3 strand (from the C-terminal V1/V2 stem): residues 201-203
data.set.13k <- data.frame (tmscore.region.beta3=data.analysis$tmscore.region.beta3)

#    (l) Beta20 strand (from the N-terminal C4 region): residues 421-423
data.set.13l <- data.frame (tmscore.region.beta20=data.analysis$tmscore.region.beta20)

#    (m) Beta21 strand (from the C-terminal C4 region): residues 433-435
data.set.13m <- data.frame (tmscore.region.beta21=data.analysis$tmscore.region.beta21)

#    (n) Items (i) through (l) combined
data.set.13n <- data.analysis[, 7088:7091]

#    (o) Items (b) through (l) combined
data.set.13o <- data.analysis[, c (7082:7087, 7093, 7088:7091)]

#    (p) All of gp120
data.set.13p <- data.frame (tmscore.region.gp120=data.analysis$tmscore.region.gp120)


# ---------------------------------------------------------------------------- #
# STEP 3:  save out for separate analysis on special vintage laptop
# ---------------------------------------------------------------------------- #


# save just the "data.set" variables, since that should be all we need
save (list=c ("data.treatment", ls ()[grep ("data.set.", ls (), fixed=T)]), file="/path/to/home/directory/results/v705_hidim_data_v1.Rdata")


# ---------------------------------------------------------------------------- #
#                                    - 30 -
# ---------------------------------------------------------------------------- #



#  1. Control variables only
#    -- Age, BMI, baseline risk score and enrollment in South Africa
#
#  2. Tier 1: Env gp120 C97ZA insert residue match/mismatch status
#    (a) V1V2-hotspot Janssen NHP challenge: sites 157-184
#    (b) IgG3 V1V2 breadth-correlates sites: sites 130, 155, 160, 161, 164, 165, 166, 170, 172, 173, 178, 179, 181, 186, 200
#    (c) RV144 V2 C1 region sites within ADCC epitope: sites 51-61, 69-78
#    (d) RV144 primary sieve sites (Rolland, Edlefsen et al. 2012): sites 169 and 181
#    (e) Signature sites for C-clade bnAbs (i.e., VRC01, 4E10, 3BNC117, and CAP256.VRC26):  sites 230, 279-282, 350, 364, 429, 432, 456, 471
#
#  3. Tier 1: Env gp70-1428 V1V2 antigen sequence residue match/mismatch status
#    -- gG3 V1V2 breadth correlates sites: sites 130, 155, 160, 161, 164, 165, 166, 170, 172, 173, 178, 179, 181, 186, 200
#
#  4. Tier 1: Env gp120 C97ZA vaccine insert-independent residue presence/absence indicators
#    (a) V1V2-hotspot Janssen NHP challenge: sites 157-184
#    (b) IgG3 V1V2 breadth-correlates sites: sites 130, 155, 160, 161, 164, 165, 166, 170, 172, 173, 178, 179, 181, 186, 200
#    (c) RV144 V2 C1 region sites within ADCC epitope: sites 51-61, 69-78
#    (d) RV144 primary sieve sites (Rolland, Edlefsen et al. 2012): sites 169 and 181
#    (e) Signature sites for C-clade bnAbs (i.e., VRC01, 4E10, 3BNC117, and CAP256.VRC26):  sites 230, 279-282, 350, 364, 429, 432, 456, 471
#
#  5. Tier 1: Presence/absence of sequons
#    (a) Glycan holes in C97ZA insert (Korber CHAVD presentation, Oct. 2021): sites 234, 392
#    (b) V2 bnAb contact glycans: sites 130, 156, 160, 234, 332, 362, 363, 392
#
#  6. Tier 1: PC-weighted Hamming distance from vaccine insert
#    (a) V1V2-hotspot Janssen NHP challenge: sites 157-184
#    (b) IgG3 V1V2 breadth correlates sites: sites 130, 155, 160, 161, 164, 165, 166, 170, 172, 173, 178, 179, 181, 186, 200
#    (c) RV144 V2 C1 region sites within ADCC epitope: sites 51-61 and 69-78
#    (d) Signature sites for C-clade bnAbs (i.e., VRC01, 4E10, 3BNC117, and CAP256.VRC26):  sites 230, 234, 279-282, 350, 364, 429, 432, 456, 471
#    (e) HVTN 505 sieve analysis (DeCamp et al.): 93 CD4bs mAb contact sites
#    (f) HVTN 505 sieve analysis (DeCamp et al.): 54 sites corresponding to 9-mers in 4 linear regions overlapping the CD4bs contact residues
#
#  7. Tier 1: PC-weighted Hamming dist from gp70-1428 V1V2 antigen
#    -- IgG3 V1V2 breadth correlates sites: sites 130, 155, 160, 161, 164, 165, 166, 170, 172, 173, 178, 179, 181, 186, 200
#
#  8. Tier 1: Other Env AA features
#    (a) Length of V1V2
#    (b) Number of sequons in V1V2
#    (c) Electrochemical charge of V2
#    (d) Cysteine count in gp160
#    (e) Indicators of single vs multiple viral lineages
#
#  9. Tier 2: Vaccine residue match/mismatch indicators in vaccine-insert region of Env, excluding Tier 1 hypothesis-driven position set
#    -- All sufficiently-variable positions in signal peptide + gp120 + N-terminal portion of gp41 included in vaccine inserts, excluding Tier 1 hypothesis-driven positions to avoid duplication
#
#  10. Tier 2: Vaccine insert-independent residue presence/absence indicators in all of Env, excluding Tier 1 hypothesis-driven position set
#    -- All sufficiently-variable positions in Env, excluding Tier 1 hypothesis-driven positions to avoid duplication (vaccine insert residues are ignored in the vaccine-insert region)
#
#  11. Tier 2: PC-weighted Hamming distance to the vaccine insert
#    (a) All of gp120
#    (b) N-terminal portion of gp41 included in vaccine inserts
#    (c) All of V5
#
#  12. Tier 2: Other AA features
#    (a) Length of V5
#    (b) Number of sequons in V5 (dichotomized as < 2 vs. >= 2)
#
#  13. Structural Features: TM score, comparing sequence with C97ZA insert
#    (a) All 15mer window TM scores
#    (b) C1 region
#    (c) V1/V2 loop
#    (d) Loop D
#    (e) C2 region
#    (f) C3 region
#    (g) C4 region
#    (h) C5 region
#    (i) Items (b) through (h) combined
#    (j) Beta2 strand (from the N-terminal V1/V2 stem): positions 121-123
#    (k) Beta3 strand (from the C-terminal V1/V2 stem): residues 201-203
#    (l) Beta20 strand (from the N-terminal C4 region): residues 421-423
#    (m) Beta21 strand (from the C-terminal C4 region): residues 433-435
#    (n) Items (i) through (m) combined
#    (o) Items (i) and (n) combined
#    (p) All of gp120


