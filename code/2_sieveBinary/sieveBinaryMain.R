# Purpose: Estimation of hazard ratio-based VE by each Tier 1 and Tier 2 binary features (in SAP Table 2); 
#          Considers primary endpoints, and the PP cohort is the analysis cohort
#          Hypothesis testing; plotting
#          Time of the first RNA positive sample used as the failure time variable
# Method:  Lunn and McNeil (1995)
# Input:   705 primary endpoints and binary features
# Output:  csv and PDF files

# refresh the workspace
rm(list=ls(all=TRUE))

# Setting directory paths -------------------------------------------------
here::i_am("README.md")
repoDir <- here::here()
dataDir <- file.path(repoDir, "data")
codeDir <- file.path(repoDir, "code/2_sieveBinary")
figureDir <- file.path(repoDir, "figures")
tableDir <- file.path(repoDir, "tables")


# Initialization
library(tidyverse)
library(plyr)
library(grid)
library(gridExtra)
library(gtable)


# load supplement
source(file.path(codeDir, "lunnMcneil.R"))
source(file.path(codeDir, "sieveBinaryUtils.R"))
source(file.path(codeDir, "forest.R"))
source(file.path(repoDir, "code/common.R"))

# load input data
sieveData <- read.csv(file.path(dataDir, datFile)) %>%
  filter(cohort=="Per-Protocol") %>%
  mutate(tx=as.numeric(armdesc=="Vaccine"), 
         lineageLabel = as.numeric(transmitted.founder.status=="multiple"),
         eventTime = hiv1fposday, 
         eventInd = hiv1event,
         numSequonV5 = as.numeric(num.sequons.v5 >= 2))

######################################################################################################
posIsAA <- colnames(sieveData)[grepl(".is.", colnames(sieveData)) & grepl("hxb2", colnames(sieveData)) & 
                                 (!grepl("is.sequon.tier1", colnames(sieveData)))]
sequonPos <- colnames(sieveData)[grepl("is.sequon.tier1", colnames(sieveData))]
posIsAAtier1PosScreenedIn <- read.csv(file.path(tableDir, paste0("posIsAAtier1posVarScreenedIn.csv")))
posIsAAtier2PosScreenedIn <- read.csv(file.path(tableDir, paste0("posIsAAtier2posVarScreenedIn.csv")))
sequonPosScreenedIn <- read.csv(file.path(tableDir, paste0("sequonPosVarScreenedIn.csv")))



# tier 1 match/mismatch
tier1binary <- colnames(sieveData)[grepl("tier1",colnames(sieveData))]
tier1posScreenedInList <-  list(posIsAA = posIsAAtier1PosScreenedIn$x[posIsAAtier1PosScreenedIn$x %in% tier1binary])
tier2posScreenedInList <- list(posIsAA = posIsAAtier2PosScreenedIn$x[!posIsAAtier2PosScreenedIn$x %in% tier1binary])
binaryMarksList <- list( tier1 = c(tier1posScreenedInList$posIsAA,
                                   sequonPosScreenedIn$x,
                                   "lineageLabel"),
                         tier2 = c(tier2posScreenedInList$posIsAA,
                                   "numSequonV5"))



# Fit competing risks survival models to obtain VE for all binary marks
binaryResultList <- list()
for(tier in c("tier1","tier2")){
  tierResultList <- list()
  for(mark in binaryMarksList[[tier]]){
    fitData <- dplyr::select(sieveData, all_of(c("eventTime", "eventInd", mark, "tx", "ind_sa")))
    # convert mark values for non-primary endpoints through tau to NA
    fitData$mark<- ifelse(fitData$eventInd==0, NA, fitData[,mark])
    # complete-case analysis, i.e., discard cases with a missing mark
    fitData <- filter(fitData, !(eventInd==1 & is.na(mark)))
    counts <- data.summary(fitData, "mark")
    #mark-specific cox regression
    fitTable <- tibble(mark = character(), inc = character(), VE = numeric(), LB = numeric(), UB = numeric(), p = numeric())
    for(i in c(1,0)){
      fitDataMS <- mutate(fitData, newStatus = as.numeric((!is.na(mark)) & mark==i))
      cox<- coxph(Surv(eventTime, newStatus) ~ tx + strata(ind_sa), data=fitDataMS)
      scox <- summary(cox)$coef
      VE <- 1 - scox[2]
      LB <- 1 - exp(scox[1] + qnorm(0.975)*scox[3])
      UB <- 1 - exp(scox[1] - qnorm(0.975)*scox[3])
      p <- scox[5]
      vaccineCases <- counts$vaccineCases[paste(i)]
      placeboCases <- counts$placeboCases[paste(i)]
      vaccineIncRate <- counts$vaccineIncRate[paste(i)]
      placeboIncRate <- counts$placeboIncRate[paste(i)]
      fitTable <- add_row(. = fitTable, 
                          mark = ifelse(i==0, "0", "1"),
                          inc = paste0(vaccineCases," (",vaccineIncRate,")", " vs. ",
                                       placeboCases," (",placeboIncRate,")"),
                          VE = VE, LB = LB, UB = UB, p = p
      )
    }
    lunnMcneilP <- with(fitData, lunnMcneilTestS(eventTime,eventInd,mark+1,tx,ind_sa))$coef[3,5]
    tierResultList[[mark]] <- list(VEtable = fitTable, diffP = lunnMcneilP)
  }
  binaryResultList[[tier]] <- tierResultList
}



# Table summary for tier 2 type 1
markResultList <- binaryResultList[["tier2"]]
WestfallYoungAdjPvalues <- read.csv(file=file.path(tableDir, paste0("WestfallYoungAdjPvalues_tier2Type1and2.csv")))
colnames(WestfallYoungAdjPvalues ) <- c("X", "p.unadj","p.FWER","p.FDR")
posScreenedInList <- tier2posScreenedInList
posIsAAVE = table.VE.residue(markResultList,marks = tier2posScreenedInList$posIsAA, WestfallYoungAdjPvalues) 
restrictedSubset <- filter(WestfallYoungAdjPvalues, p.unadj <= 0.05)$X
pvalues <- ldply(markResultList, function(x)x$diffP)
colnames(pvalues) <- c("mark", "unadjustedP")
pvaluesSubset <- filter(pvalues, unadjustedP <= 0.05)
write.csv(posIsAAVE, file.path(tableDir,"sieveBinary", paste0("tier2posIsAAVE.csv")), row.names = FALSE)


# Table summary for tier 1 type 1 and 3
markResultList <- binaryResultList[["tier1"]]
WestfallYoungAdjPvalues <- read.csv(file=file.path(tableDir, paste0("WestfallYoungAdjPvalues_tier1Type1to4.csv")))
colnames(WestfallYoungAdjPvalues) <- c("X", "p.unadj", "p.FWER", "p.FDR")

posScreenedInList <- tier1posScreenedInList
posIsAAVE = table.VE.residue(markResultList, marks = tier1posScreenedInList$posIsAA, WestfallYoungAdjPvalues)
tier1Positions <- colnames(sieveData)[grepl("tier1",colnames(sieveData))]
tier1aaPos <- laply(tier1Positions[1:68],function(list)strsplit(list,split="\\.")[[1]][[2]])
subTypeA <- tier1aaPos[loc(tier1aaPos,"157") : loc(tier1aaPos,"184")]
subTypeB <- c("130", "155","160","161","164","165","166","170","172","173","178","179","181","186","200")
subTypeC <- tier1aaPos[c(loc(tier1aaPos,"51") : loc(tier1aaPos,"61"), loc(tier1aaPos,"69") : loc(tier1aaPos,"78"))]
subTypeD <- c("169", "181")
subTypeE <- c("230", tier1aaPos[c(loc(tier1aaPos,"279") : loc(tier1aaPos,"282"))], "350","364","429","432","456","471")

positionExtract <- function(positions)laply(positions, function(list)strsplit(list, split="\\$")[[1]][[1]])

subTypeList <- list("subTypeA" = subTypeA,"subTypeB" = subTypeB,"subTypeC" = subTypeC,
                    "subTypeD" = subTypeD,"subTypeE" = subTypeE)
for(subType in c("subTypeA", "subTypeB", "subTypeC", "subTypeD", "subTypeE")){
  posIsAAVEsubType <- subset(posIsAAVE, positionExtract(position) %in% subTypeList[[subType]])
  write.csv(posIsAAVEsubType, file.path(tableDir,"sieveBinary", paste0("tier1posIsAAVE_", subType,".csv")),row.names = FALSE)
}



# forestplot for tier 1 multiplicity of lineage
fit <- binaryResultList$tier1$lineageLabel
VEtable <- fit$VEtable
VE <- tibble(feature = character(), cases = character(), VE = character(),mean = numeric(), ll = numeric(), ul = numeric(),
             P = character(), diffP = character())
VE <- add_row(VE, feature = "", cases = "", VE = "" , mean = NA, ll = NA, ul = NA, P = "",
              diffP = format.p(fit$diffP))
VE <- add_row(VE, feature  = "Multiple Lineages", cases = VEtable$inc[1],
              VE = paste0(format.VE(VEtable$VE[1]), " (", format.VE(VEtable$LB[1]),", ",format.VE(VEtable$UB[1]),")") ,
              mean = VEtable$VE[1]*100, ll = VEtable$LB[1]*100, ul = VEtable$UB[1]*100,
              P = format.p(VEtable$p[1]), diffP = "")
VE <- add_row(VE, feature = "Single Lineage", cases = VEtable$inc[2],
              VE = paste0(format.VE(VEtable$VE[2]), " (", format.VE(VEtable$LB[2]),", ",format.VE(VEtable$UB[2]),")") ,
              mean = VEtable$VE[2]*100, ll = VEtable$LB[2]*100, ul = VEtable$UB[2]*100,
              P = format.p(VEtable$p[2]), diffP = "")
header.plot <- c("Lineage\nMultiplicity", "No. of Cases (V vs. P)\n(Incidence per 100 PYRs)", "VE (%) (95% CI)", "Two-sided\nP-value", "Two-sided Differential VE\nUnadjusted P-value")
p <- forestWrap ( VE, c(4,5,6), header.plot, xlim = c(-100, 100), 
                  x_ticks_at = c(-100, -50, 0, 50, 100), xlab = "VE (%) (95% CI)", 
                  ci_cols = rep(c("black","royalblue", "darkred"),dim(VE)[1]/3),
                  nrows_underline = NULL,
                  insertHeaderText = NULL)

ggplot2::ggsave(filename = "VE_singleVsMultipleFounder.pdf", 
                plot = p, 
                path = figureDir,
                dpi = 320,
                width = 9, height = 2,units = "in") 


# forestplot for residue scanning at position 364
WestfallYoungAdjPvalues <- read.csv(file=file.path(tableDir, paste0("WestfallYoungAdjPvalues_tier1Type1to4.csv")))
colnames(WestfallYoungAdjPvalues) <- c("X", "p.unadj","p.FWER","p.FDR")
marks <- c("hxb2.364.is.S.tier1","hxb2.364.is.P.tier1")
VE <- tibble(mark = character(), position = character(),residue = character(), cases = character(), VE = character(), mean = numeric(),
             ll = numeric(), ul = numeric(), P = character(), 
             diffP = character(), FWER = character(), Qvalue = character() )
for(mark in marks){
  VEtable <- binaryResultList$tier1[[mark]]$VEtable
  y = strsplit(mark,split="\\.")[[1]]
  pos <- y[2]
  residue <- y[4]
  unadjp <- binaryResultList$tier1[[mark]]$diffP
  fdr = subset(WestfallYoungAdjPvalues, X==mark)$"p.FDR"
  holm = subset(WestfallYoungAdjPvalues, X==mark)$"p.FWER"
  Qvalue = format.p(fdr)
  if(fdr <= 0.2 & unadjp <= 0.05){Qvalue = paste0(Qvalue,"*")}
  VE <- add_row(VE, 
                mark = mark,
                position = pos,
                residue = "",
                cases = "",
                VE = "",
                mean = NA, ll = NA, ul = NA, P = "",
                diffP = format.p(binaryResultList$tier1[[mark]]$diffP), 
                FWER = format.p(holm), Qvalue = Qvalue)
  VE <- add_row(VE, 
                mark = mark,
                position = "",
                residue = residue,
                cases = VEtable$inc[1],
                VE = paste0(format.VE(VEtable$VE[1]), " (", format.VE(VEtable$LB[1]),", ",format.VE(VEtable$UB[1]),")") ,
                mean = VEtable$VE[1]*100, ll = VEtable$LB[1]*100, ul = VEtable$UB[1]*100,P = format.p(VEtable$p[1]), 
                diffP = "", 
                FWER = "",Qvalue = "")
  VE <- add_row(VE, 
                mark = mark,
                position = "",
                residue = paste0 ("Not ", residue),
                cases = VEtable$inc[2],
                VE = paste0(format.VE(VEtable$VE[2]), " (", format.VE(VEtable$LB[2]),", ",format.VE(VEtable$UB[2]),")") ,
                mean = VEtable$VE[2]*100, ll = VEtable$LB[2]*100, ul = VEtable$UB[2]*100,P = format.p(VEtable$p[2]), 
                diffP = "", 
                FWER = "", Qvalue = "")
}

header.plot <- c("Env 364\nResidue", "No. of Cases (V vs. P)\n(Incidence per 100 PYRs)", "VE (%) (95% CI)", "Two-sided\nP-value", "P-value", "FWER\nP-value","Q-value")
p <- forestWrap ( VE[,-c(1,2)], c(4,5,6), header.plot, xlim = c(-100, 100), 
                  x_ticks_at = c(-100, -50, 0, 50, 100), xlab = "VE (%) (95% CI)", 
                  ci_cols = rep(c("black","royalblue", "darkred"),dim(VE)[1]/3),
                  nrows_underline = NULL,
                  insertHeaderText = "Two-sided Differential VE",
                  insertHeaderTextCol = c(6:8))

ggplot2::ggsave(filename = paste0("tier1Type1to4site364.pdf"), 
                plot = p, 
                path = figureDir,
                dpi = 320,
                width = 10, height = 0.25*(length(VE$mean)+4),units = "in") 


# forestplot for the presence of sequon at position 362 and 363
sequonsVE <- table.VE.sequon.wo.adj(binaryResultList$tier1, marks =  c("hxb2.362.is.sequon.tier1", "hxb2.363.is.sequon.tier1"), forestplot = TRUE)
header.plot <- c("PNGS\nStarting\nPosition", "PNGS\nPresence", "No. of Cases (V vs. P)\n(Incidence per 100 PYRs)", "VE (%) (95% CI)", "Two-sided\nP-value", "Differential VE\nUnadjusted P-value")
p <- forestWrap (sequonsVE[,-c(1)], c(5,6,7), header.plot, xlim = c(-100, 100), 
                  x_ticks_at = c(-100, -50, 0, 50, 100), xlab = "VE (%) (95% CI)", 
                  ci_cols = rep(c("black", "royalblue", "darkred"), dim(sequonsVE)[1]/3),
                  nrows_underline = NULL)

ggplot2::ggsave(filename = paste0("posthoc_sequon362_363.pdf"), 
                plot = p, 
                path = figureDir,
                dpi = 320,
                width = 10, height = 0.25*(length(sequonsVE$mean)+4),units = "in") 


# forestplot for the presence of sequon at positions other than 362 and 363
markResultList <- binaryResultList[["tier1"]]
WestfallYoungAdjPvalues <- read.csv(file=file.path(tableDir, paste0("WestfallYoungAdjPvalues_tier1Type1to4.csv")))
colnames(WestfallYoungAdjPvalues ) <- c("X", "p.unadj", "p.FWER", "p.FDR")
sequonPosScreenedIn2 <- sequonPosScreenedIn$x[!sequonPosScreenedIn$x %in% c("hxb2.362.is.sequon.tier1", "hxb2.363.is.sequon.tier1")]
sequonsVE_forestTable <- table.VE.sequon(markResultList, marks = sequonPosScreenedIn2, WestfallYoungAdjPvalues, forestplot = TRUE)
sequonsVE_forestTable$position[sequonsVE_forestTable$position == "234"] = "234*"
sequonsVE_forestTable$position[sequonsVE_forestTable$position == "392"] = "392*"
header.plot <- c("PNGS\nStarting\nPosition", "PNGS\nPresence", "No. of Cases (V vs. P)\n(Incidence per 100 PYRs)", "VE (%) (95% CI)", "Two-sided\nP-value", "P-value", "FWER\nP-value","Q-value")
p <- forestWrap ( sequonsVE_forestTable[,-1], c(5,6,7), header.plot, xlim = c(-30, 110), 
                  x_ticks_at = c(-20, 0, 20, 40, 60, 80), xlab = "VE (%) (95% CI)", 
                  ci_cols = rep(c("black","royalblue", "darkred"), dim(sequonsVE_forestTable)[1]/3),
                  nrows_underline = NULL,
                  insertHeaderText = "Two-sided Differential VE",
                  insertHeaderTextCol = c(7:9),
                  footnote = "\n* Glycan hole in the C97ZA insert sequence")

ggplot2::ggsave(filename = paste0("tier1sequon.pdf"), 
                plot = p, 
                path = figureDir,
                dpi = 320,
                width = 10, height = 0.25*(length(sequonsVE_forestTable$mean)+7),units = "in") 



# Forest plot summary for tier 2 type 3: number of sequons in V5
WestfallYoungAdjPvalues <- read.csv(file=file.path(tableDir, paste0("WestfallYoungAdjPvalues_tier2Type3and4.csv")))
colnames(WestfallYoungAdjPvalues ) <- c("X", "p.unadj","p.FWER","p.FDR")
markResultList <- binaryResultList[["tier2"]]
VEtable <- markResultList[["numSequonV5"]]$VEtable
adjustedP <- subset(WestfallYoungAdjPvalues, X == "numSequonV5")
numSequonV5VE_forestplot <- tibble(mark = character(), cases = character(), VE = character(),mean = numeric(), ll = numeric(), ul = numeric(),
                                   P = character(), diffP = character(), FWER = character(), Qvalue = character())
numSequonV5VE_forestplot <- add_row(numSequonV5VE_forestplot, mark = "", cases = "", VE = "", mean = NA, ll = NA, ul = NA,
                                    P = "", diffP = format.p( adjustedP$p.unadj), FWER = format.p(adjustedP$p.FWER), Qvalue = format.p(adjustedP$p.FDR))
numSequonV5VE_forestplot <- add_row(numSequonV5VE_forestplot, mark = "<= 1", cases = VEtable$inc[2], VE = numSequonV5VE$VE0, mean = VEtable$VE[2]*100, ll = VEtable$LB[2]*100, ul = VEtable$UB[2]*100,
                                    P = format.p(VEtable$p[2]), diffP = "", FWER = "", Qvalue = "")
numSequonV5VE_forestplot <- add_row(numSequonV5VE_forestplot, mark = ">= 2", cases = VEtable$inc[1], VE = numSequonV5VE$VE1, mean = VEtable$VE[1]*100, ll = VEtable$LB[1]*100, ul = VEtable$UB[1]*100,
                                    P = format.p(VEtable$p[1]), diffP = "", FWER = "", Qvalue = "")
header.plot <- c("No. of\nPNGS motifs\nin V5", "No. of Cases (V vs. P)\n(Incidence per 100 PYRs)", "VE (%) (95% CI)", "Two-sided\nP-value", "P-value", "FWER\nP-value","Q-value")
p <- forestWrap ( numSequonV5VE_forestplot, c(4,5,6), header.plot, xlim = c(-30, 100), 
                  x_ticks_at = c(-20, 0, 20, 40, 60, 80), xlab = "VE (%) (95% CI)", 
                  ci_cols = rep(c("black","royalblue", "darkred"),dim(numSequonV5VE_forestplot)[1]/3),
                  nrows_underline = NULL,
                  insertHeaderText = "Two-sided Differential VE",
                  insertHeaderTextCol = c(6:8))

ggplot2::ggsave(filename = "numSequonV5VE_forestplot.pdf", 
                plot = p, 
                path = figureDir,
                dpi = 320,
                width = 10, height = 0.25*(length(numSequonV5VE_forestplot$mean)+5), units = "in")  




