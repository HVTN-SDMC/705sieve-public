# Purpose: Westfall and Young permutation-based multiplicity adjustment for sieve test p-values
#          The SAP states that p-values are calculated for the vaccine vs. placebo comparison.
# Method:  Westfall and Young (1993)
#          Juraska and Gilbert (Biometrics, 2013)
#          R package sievePH, version 1.0.3 on CRAN
#          Lunn and McNeil (1995)
# Input:   705 primary endpoints,  sequence features specified in the sieve SAP that are to be adjusted in multiple comparison
# Output:  A list with each component being a vector of sieve test p-values for the analyzed marks from a single permutation of the mark variable


# refresh the workspace
rm(list=ls(all=TRUE))

# Setting directory paths -------------------------------------------------
dataDir <- "/file/path/data"
codeDir <- "/file/path/code/westfallYoung"
figureDir <- "/file/path/figures/westfallYoung"
tableDir <- "/file/path/tables/westfallYoung"
preScreentableDir <- "/file/path/tables/preScreen"

# Initialization
library(tidyverse)
library(plyr)
library(sievePH)
# load supplement
source(file.path(codeDir, "lunnMcneil.R"))
source(file.path(codeDir, "p.adj.perm2.R"))
source(file.path(codeDir, "common.R"))

# load input data
sieveData <- read.csv(file.path(dataDir, datFile)) %>%
  filter(cohort=="Per-Protocol") %>%
  mutate(tx=as.numeric(armdesc=="Vaccine"), 
         lineageLabel = as.numeric(transmitted.founder.status=="multiple"),
         eventTime = hiv1fposday, 
         eventInd = hiv1event,
         numSequonV5 = as.numeric(num.sequons.v5 >= 2))

######################################################################################################
posIsAA <- colnames(sieveData)[grepl(".is.",colnames(sieveData)) & grepl("hxb2",colnames(sieveData)) & (!grepl("is.sequon.tier1",colnames(sieveData)))]
sequonPos <- colnames(sieveData)[grepl("is.sequon.tier1",colnames(sieveData))]
hdist.zspace <- colnames(sieveData)[grepl("hdist.zspace",colnames(sieveData))]
hdist_ab_no364 <- grep("no364", hdist.zspace, value = TRUE)
hdist.zspace <- hdist.zspace[!hdist.zspace%in% hdist_ab_no364 ]

# load binary marks that pass the screening filter
posIsAAtier1PosScreenedIn <- read.csv(file.path(preScreentableDir , paste0("posIsAAtier1posVarScreenedIn.csv")))
posIsAAtier2PosScreenedIn <- read.csv(file.path(preScreentableDir , paste0("posIsAAtier1posVarScreenedIn.csv")))

sequonPosScreenedIn <- read.csv(file.path(preScreentableDir , paste0("sequonPosVarScreenedIn.csv")))
#Exclude sequon at postions 362 and 363 since they were sensitivity analysis
sequonPosScreenedIn2 <- sequonPosScreenedIn$x[!sequonPosScreenedIn$x %in% c("hxb2.362.is.sequon.tier1", "hxb2.363.is.sequon.tier1")]


# tier 1 binary
tier1binary <- colnames(sieveData)[grepl("tier1",colnames(sieveData))]

# specify the set of continuous marks in tier 1 and tier 2; not all marks are adjusted for P-values
continousMarksList <- list(tier1 = c("hdist.zspace.mos1.v2","hdist.zspace.mos2.v2","hdist.zspace.c97za.v2" ,
                                     "hdist.zspace.mos1.v2_ab","hdist.zspace.mos2.v2_ab","hdist.zspace.c97za.v2_ab",
                                     "hdist.zspace.mos1.adcc","hdist.zspace.mos2.adcc","hdist.zspace.c97za.adcc",
                                     "hdist.zspace.mos1.c_ab","hdist.zspace.mos2.c_ab","hdist.zspace.c97za.c_ab",
                                     "hdist.zspace.c97za.hvtn505.cd4bs.antibody", "hdist.zspace.c97za.hvtn505.cd4bs.kmer" ,
                                     "hdist.zspace.1428v1v2",
                                     "length.v1v2",
                                     "num.sequons.v1v2",
                                     "charge.v2",
                                     "cysteine.count"),
                           tier2 = c("hdist.zspace.mos1.gp120","hdist.zspace.mos2.gp120","hdist.zspace.c97za.gp120",
                                     "hdist.zspace.mos1.gp41","hdist.zspace.mos2.gp41","hdist.zspace.c97za.gp41",
                                     "hdist.zspace.mos1.v5","hdist.zspace.mos2.v5","hdist.zspace.c97za.v5",
                                     "length.v5"                                        )
)


# multiple testing adjust groups
adjustList <- list(tier2Type1and2 = c(posIsAAtier2PosScreenedIn$x[!posIsAAtier2PosScreenedIn$x %in% tier1binary]),
                   tier2Type3and4 = c("hdist.zspace.c97za.gp120", "hdist.zspace.c97za.v5", "hdist.zspace.c97za.gp41", "numSequonV5", "length.v5"),
                   tier1Type1to4 = c(posIsAAtier1PosScreenedIn$x[posIsAAtier1PosScreenedIn$x% in% tier1binary], sequonPosScreenedIn2),
                   tier1Type5to7= c("hdist.zspace.c97za.v2", "hdist.zspace.c97za.v2_ab","hdist.zspace.c97za.adcc","hdist.zspace.c97za.c_ab",
                                    "hdist.zspace.c97za.hvtn505.cd4bs.antibody", "hdist.zspace.c97za.hvtn505.cd4bs.kmer" ,
                                    "hdist.zspace.1428v1v2",
                                    "length.v1v2","num.sequons.v1v2","charge.v2","cysteine.count")
)



# Compute p-values from data sets with resampled marks --------------------
nPerm <- 1000
for(adjustGroup in c("tier1Type1to4", "tier1Type5to7", "tier2Type1and2", "tier2Type3and4")){
  marks <- adjustList[[adjustGroup]]
  sieveDataSub <- select(sieveData, all_of(c("eventTime", "eventInd", "tx","ind_sa", marks)))
  # get the p-values for individual permutations of the observed marks
  # the first vector in 'pvals' stores unadjusted p-values based on original data
  pvals <- lapply(1:(nPerm + 1), function(seed){
    set.seed(seed)
    
    # permute the marks observed or missing in cases
    idx <- sample(1:sum(sieveDataSub$eventInd))
    
    # 'pvals1' is a vector of p-values (one for each mark variable) for the single permutation
    pvals1 <- sapply(1:length(marks), function(i){
      data1 <- select(sieveDataSub,all_of(c("eventTime", "eventInd", "tx", "ind_sa", marks[i])))
      colnames(data1) <- c("eventTime", "eventInd", "tx", "ind_sa", "mark")
      
      # convert mark values for non-primary endpoints to NA
      data1$mark <- ifelse(data1$eventInd == 0, NA, data1$mark)
      
      # apply the permutation
      if (seed > 1){
        data1$mark[data1$eventInd == 1] <- data1$mark[data1$eventInd == 1][idx]  
      }
      
      # complete-case analysis, i.e., discard cases with a missing mark
      data1 <- subset(data1, !(eventInd==1 & is.na(mark)))
      
      # if the mark is quantitative
      if (marks[i] %in% unlist(continousMarksList)){
        # fit the mark-specific HR model
        markRng <- range(data1$mark, na.rm=TRUE)
        markGrid <- seq(markRng[1], markRng[2], length.out=200)
        fit <- with(data1, sievePH(eventTime=eventTime, eventInd=eventInd, mark=mark, tx=tx, strata=ind_sa))
        if(marks[i] %in% hdist.zspace){
          sfit <- summary(fit, markGrid=markGrid, sieveAlternative="oneSided")
          # double 1-sided Wald test of {H0: PE(v) constant for all v}
          return(min(2*sfit$pWald.HRconstant.1sided,1))
        }else{
          sfit <- summary(fit, markGrid=markGrid, sieveAlternative="twoSided")
          # 2-sided Wald test of {H0: PE(v) constant for all v}
          return(sfit$pWald.HRconstant.2sided) 
        }
      } else {
        sfit <- with(data1, lunnMcneilTestS(eventTime, eventInd, mark + 1, tx, stratVar=ind_sa))
        return(sfit$coef[3, 5])
      }
    })
    
    names(pvals1) <- marks
    return(pvals1)
  })
  
  pvals <- do.call(rbind, pvals)
  # Apply Westfall and Young (1993) to obtain adjusted p-values -------------
  pvals.adj <- p.adj.perm2(p.unadj=pvals[1, ], p.perms=pvals[-1, ], alpha=1)
  mark <- as.vector(row.names(pvals.adj))
  npvals.adj <- data.frame(mark,pvals.adj )
  write.csv(npvals.adj, file=file.path(tableDir, paste0("WestfallYoungAdjPvalues_",adjustGroup,".csv")), row.names=FALSE)
  save(pvals, file=file.path(tableDir, paste0("WestfallYoungPermPvalues_",adjustGroup,".RData")))
}




