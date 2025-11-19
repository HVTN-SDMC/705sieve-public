

# ---------------------------------------------------------------------------- #
# STEP 0:  set up our working environment
# ---------------------------------------------------------------------------- #


# refresh workspace
rm (list=ls ())
set.seed (378)

# load required packages
library ("diproperm")


# ---------------------------------------------------------------------------- #
# STEP 1:  load and process our indata
# ---------------------------------------------------------------------------- #


# where does our data live?
path.home <- "/path/to/home/directory"

# load Rdata object
load (file.path (path.home, "v705_hidim_data_v1.Rdata"))

# we need a function to rescale our datasets
rescale <- function (data.in) {

  # process one column at a time
  for (col.num in 1:ncol (data.in)) {
    if (var (data.in[, col.num], na.rm=T) == 0) {
      data.in[, col.num] <- rep (1, nrow (data.in))
    } else {
      data.tmp <- data.in[, col.num]
      data.in[!is.na (data.in[, col.num]), col.num] <- (data.in[!is.na (data.in[, col.num]), col.num] - min (data.in[!is.na (data.in[, col.num]), col.num])) / max (data.in[!is.na (data.in[, col.num]), col.num] - min (data.in[!is.na (data.in[, col.num]), col.num]))
    }
  }

  # done, give it up
  return (data.in)
}

# we need a list of our set names/numbers
analysis.sets <- ls ()[grep ("data.set", ls (), fixed=T)]


# ---------------------------------------------------------------------------- #
# STEP 2:  get everything ready
# ---------------------------------------------------------------------------- #


# our class is always going to be the treatment assignment; format as required
class <- data.treatment
class[class == 0] <- -1

# fix the categorical data in set 8e
data.set.8e[is.na (data.set.8e[, 1]), 1] <- "0"
data.set.8e[data.set.8e[, 1] == "single", 1] <- "1"
data.set.8e[data.set.8e[, 1] == "multiple", 1] <- "2"
data.set.8e[, 1] <- as.numeric (data.set.8e[, 1]) - 1
data.set.8e[data.set.8e[, 1] == -1, 1] <- NA

# rescale all of our datasets
for (set in analysis.sets) {
  assign (set, rescale (eval (as.name (paste (set)))))
}

# we need to create data supersets, to combine subsets into a whole, when 
# appropriate
data.set.2 <- cbind (data.set.2a, data.set.2b, data.set.2c, data.set.2d, data.set.2e)
data.set.4 <- cbind (data.set.4a, data.set.4b, data.set.4c, data.set.4d, data.set.4e)
data.set.5 <- cbind (data.set.5a, data.set.5b)
data.set.6 <- cbind (data.set.6a, data.set.6b, data.set.6c, data.set.6d, data.set.6e, data.set.6f)
data.set.8 <- cbind (data.set.8a, data.set.8b, data.set.8c, data.set.8d, data.set.8e)
data.set.11 <- cbind (data.set.11a, data.set.11b, data.set.11c)
data.set.12 <- cbind (data.set.12a, data.set.12b)

# get a new list of our analysis set names/numbers
analysis.sets <- ls ()[grep ("data.set", ls (), fixed=T)]
analysis.sets.hybrid <- analysis.sets[grep ("13", analysis.sets) * -1][-1]

# set up a data frame for results
results <- data.frame (set=c ("null.set", analysis.sets, paste (analysis.sets.hybrid, "hybrid", sep="."), "data.set.all"), p.val=NA, num.vars=NA, matrix (NA, nrow=length (grep ("data.set", ls (), fixed=T)) + length (analysis.sets.hybrid) + 2, ncol=20))
num.loading.vars <- 10
names (results)[4:(3 + (num.loading.vars * 2))] <- as.vector (sapply (1: num.loading.vars, function (x) return (c (paste0 ("loading.", x, ".var"), paste0 ("loading.", x, ".val")))))


# ---------------------------------------------------------------------------- #
# STEP 3a:  run analysis on null set
# ---------------------------------------------------------------------------- #


# generate our null data
x <- matrix (0, nrow=length (class), ncol=1, byrow=T)
x[class == 1] <- c (rep (1, floor (sum (class == 1) / 2)), rep (2, ceiling (sum (class == 1) / 2)))
x[class == -1] <- c (rep (1, floor (sum (class == -1) / 2)), rep (2, ceiling (sum (class == -1) / 2)))
x <- cbind (x, x)
colnames (x) <- c ("null.1", "null.2")
y <- class

# run analysis, save p-value, number of variables and top ten loadings
dpp.tmp <- DiProPerm (x, y, B=1000, classifier="md", univ.stat="md")
results[results$set == "null.set", "p.val"] <- dpp.tmp$pvalue
results[results$set == "null.set", "num.vars"] <- ncol (x)
loadings.tmp <- loadings (dpp.tmp, loadnum = 10)
loadings.num <- 10 - sum (is.na (loadings.tmp$index))
results[results$set == "null.set", 3 + seq (1, (loadings.num * 2) - 1, 2)] <- colnames (x)[loadings.tmp$index][1:loadings.num]
results[results$set == "null.set", 4 + seq (1, (loadings.num * 2) - 1, 2)] <- loadings.tmp[1:loadings.num, 2]

# save plot
pdf (file="/path/to/home/directory/results/data.set.null.pdf", height=8, width=6)
plotdpp (dpp.tmp)
dev.off ()


# ---------------------------------------------------------------------------- #
# STEP 3b:  run our analysis of individual sets
# ---------------------------------------------------------------------------- #


# loop through each data set (minus the "all features" analysis)
for (set.tmp in analysis.sets) {

  # create a temporary verison of our response variable
  y <- class

  # initialize feature block
  # when running set.1 (background variables), run them alone; otherwise add
  # the background variables to the set of interest
  if (set.tmp == "data.set.1") {
    x <- as.matrix (get (set.tmp))
  } else {
    x <- as.matrix (cbind (data.set.1, get (set.tmp)))
    y <- y[rowSums (is.na (x)) == 0]
    x <- na.omit (x)
  }

  # run analysis, save p-value, number of variables and top ten loadings
  dpp.tmp <- DiProPerm (x, y, B=1000, classifier="md", univ.stat="md")
  results[results$set == set.tmp, "p.val"] <- dpp.tmp$pvalue
  results[results$set == set.tmp, "num.vars"] <- ncol (x)
  loadings.tmp <- loadings (dpp.tmp, loadnum = 10)
  loadings.num <- 10 - sum (is.na (loadings.tmp$index))
  results[results$set == set.tmp, 3 + seq (1, (loadings.num * 2) - 1, 2)] <- colnames (x)[loadings.tmp$index][1:loadings.num]
  results[results$set == set.tmp, 4 + seq (1, (loadings.num * 2) - 1, 2)] <- loadings.tmp[1:loadings.num, 2]

  # save plot
  pdf (file=paste0 ("/path/to/home/directory/results/", set.tmp, ".pdf"), height=8, width=6)
  plotdpp (dpp.tmp)
  dev.off ()
}


# ---------------------------------------------------------------------------- #
# STEP 3c:  run analyses with sequence/structural hybrid sets
# ---------------------------------------------------------------------------- #


# loop through each data set (minus the "all features" analysis)
for (set.tmp in analysis.sets.hybrid) {

  # create a temporary verison of our response variable
  y <- class

  # initialize feature block
  # when running set.1 (background variables), run them alone; otherwise add
  # the background variables to the set of interest
  x <- as.matrix (cbind (data.set.1, get (set.tmp), data.set.13o))
  y <- y[rowSums (is.na (x)) == 0]
  x <- na.omit (x)

  # run analysis, save p-value, number of variables and top ten loadings
  dpp.tmp <- DiProPerm (x, y, B=1000, classifier="md", univ.stat="md")
  results[results$set == paste (set.tmp, "hybrid", sep="."), "p.val"] <- dpp.tmp$pvalue
  results[results$set == paste (set.tmp, "hybrid", sep="."), "num.vars"] <- ncol (x)
  loadings.tmp <- loadings (dpp.tmp, loadnum = 10)
  loadings.num <- 10 - sum (is.na (loadings.tmp$index))
  results[results$set == paste (set.tmp, "hybrid", sep="."), 3 + seq (1, (loadings.num * 2) - 1, 2)] <- colnames (x)[loadings.tmp$index][1:loadings.num]
  results[results$set == paste (set.tmp, "hybrid", sep="."), 4 + seq (1, (loadings.num * 2) - 1, 2)] <- loadings.tmp[1:loadings.num, 2]

  # save plot
  pdf (file=paste0 ("/path/to/home/directory/results/", set.tmp, ".hybrid.pdf"), height=8, width=6)
  plotdpp (dpp.tmp)
  dev.off ()
}


# ---------------------------------------------------------------------------- #
# STEP 3d:  run our analysis on all features
# ---------------------------------------------------------------------------- #


# set up data objects
data.set.all <- cbind (data.set.1, data.set.2a, data.set.2b, data.set.2c, data.set.2d,
                       data.set.2e, data.set.3, data.set.4a, data.set.4b, data.set.4c,
                       data.set.4d, data.set.4e, data.set.5a, data.set.5b, data.set.6a,
                       data.set.6b, data.set.6c, data.set.6d, data.set.6e, data.set.6f,
                       data.set.7, data.set.8a, data.set.8b, data.set.8c, data.set.8d,
                       data.set.8e, data.set.9, data.set.10, data.set.11a, data.set.11b,
                       data.set.11c, data.set.12a, data.set.12b, data.set.13a, data.set.13o, 
                       data.set.13p)
x <- as.matrix (data.set.all)
y <- class
y <- y[rowSums (is.na (x)) == 0]
x <- na.omit (x)

# run analysis, save p-value and loadings
dpp.tmp <- DiProPerm (x, y, B=1000, classifier="md", univ.stat="md")
results[results$set == "data.set.all", "p.val"] <- dpp.tmp$pvalue
results[results$set == "data.set.all", "num.vars"] <- ncol (x)
loadings.tmp <- loadings (dpp.tmp, loadnum = 10)
loadings.num <- 10 - sum (is.na (loadings.tmp$index))
results[results$set == "data.set.all", 3 + seq (1, (loadings.num * 2) - 1, 2)] <- colnames (x)[loadings.tmp$index][1:loadings.num]
results[results$set == "data.set.all", 4 + seq (1, (loadings.num * 2) - 1, 2)] <- loadings.tmp[1:loadings.num, 2]

# save plot
pdf (file="/path/to/home/directory/results/data.set.all.pdf", height=8, width=6)
plotdpp (dpp.tmp)
dev.off ()


# ---------------------------------------------------------------------------- #
# STEP 4:  export our results
# ---------------------------------------------------------------------------- #


# save results into outfile
write.csv (results, file="/path/to/home/directory/results/v705_dpp_results_v2b.csv", row.names=F)

# immortalize our workspace
save.image (file="/path/to/home/directory/results/v705_dpp_v2b_image.Rdata")


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
#    (n) Items (j) through (m) combined
#    (o) Items (i) and (n) combined
#    (p) All of gp120


