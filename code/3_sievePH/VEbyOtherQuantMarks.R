# Purpose: Estimation of hazard ratio-based VE by Tier 1 feature type 4 and Tier 2 feature type 3 
#          (except number of sequons in V5, which is dichotomized)
#          Considers primary endpoints, and the PP cohort is the analysis cohort
#          Hypothesis testing; plotting
#          Time of the first RNA positive sample used as the failure time variable
# Method:  Juraska andn Gilbert (Biometrics, 2013)
#          R package sievePH, version 1.0.3 on Github
# Input:   705 primary endpoints and mark variables specified above
# Output:  PDF files, each containing a single plot

# refresh the workspace
rm(list=ls(all=TRUE))

here::i_am("README.md")
repoDir <- here::here()
dataDir <- file.path(repoDir, "data")
codeDir <- file.path(repoDir, "code/3_sievePH")
figureDir <- file.path(repoDir, "figures")
tableDir <- file.path(repoDir, "tables")

# Initialization
library(tidyverse)
library(sievePH) 

# load supplement
source(file.path(codeDir, "ggplot.summary.sievePH.R"))
source(file.path(repoDir, "code/common.R"))

# load input data
sieveData  <- read_csv(file.path(datDir, datFile)) %>%
  filter(cohort=="Per-Protocol") %>%
  select(armdesc, eventTime=hiv1fposday, eventInd=hiv1event, ind_sa, starts_with("length"), num.sequons.v1v2, 
         charge.v2, cysteine.count) %>%
  mutate(tx=as.numeric(armdesc=="Vaccine"))

######################################################################################################
# mark variable names
markName <- colnames(sieveData )[-c(1:4, ncol(sieveData))]

# axis labels
xLab <- case_when(markName=="length.v1v2" ~ "Length of V1V2",
                  markName=="length.v5" ~ "Length of V5 Loop",
                  markName=="num.sequons.v1v2" ~ "Number of PNGS Motifs in V1V2",
                  markName=="charge.v2" ~ "Electrochemical Charge of V2",
                  markName == "cysteine.count" ~ "Cysteine Count in gp160")

# strings included in file names
markFileString <- gsub(".", "_", markName, fixed=TRUE)

# initialize vectors of p-values
pHRconstancy <- pHRunity <- NULL

# for each Hamming distance
for (i in 1:length(markName)){
  
  dat1 <- select(sieveData, tx, eventTime, eventInd, mark=markName[i], ind_sa)
  
  # convert mark values for non-primary endpoints to NA
  dat1$mark <- ifelse(dat1$eventInd==0, NA, dat1$mark)
  
  # complete-case analysis, i.e., discard cases with a missing mark
  dat1 <- filter(dat1, !(eventInd==1 & is.na(mark)))
  
  # fit the mark-specific HR model
  if (markName[i] %in% c("length.v1v2")){  # plot a continuous VE curve
    markRng <- range(dat1$mark, na.rm=TRUE)
    markGrid <- seq(markRng[1], markRng[2], length.out=500)  
  } else {  # plot VE estimates as discrete points
    markGrid <- sort(unique(na.omit(dat1$mark)))
  }
  fit <- with(dat1, sievePH(eventTime=eventTime, eventInd=eventInd, mark=mark, tx=tx, strata=ind_sa))
  sfit <- summary(fit, markGrid=markGrid, sieveAlternative="twoSided")
  
  # 2-sided Wald test of {H0: VE(v)=0 for all v}
  pHRunity <- c(pHRunity, sfit$pWald.HRunity.2sided)
  
  # 2-sided Wald test of {H0: VE(v)=VE for all v} against {H1: VE(v) varies with v}
  pHRconstancy <- c(pHRconstancy, sfit$pWald.HRconstant.2sided)
  
  pdf(file.path(figDir, paste0("705_sievePH_VE_", markFileString[i], "_PP.pdf")), width=0.9*7, height=0.9*6.3)
  xLim <- range(dat1$mark, na.rm=TRUE)
  
  plotHeights <- c(0.32, 0.68)
  ggsave.width <- 0.7 * 7.3
  ggsave.height <- 0.7 * 6.5
  title <- NULL
  pval.filepath <- file.path(tableDir, paste0("WestfallYoungAdjPvalues_tier1Type5to7.csv"))
  p.df <- read_csv(pval.filepath, show_col_types=FALSE)
  p.df1 <- filter(p.df, mark%in%markName)
  
  pval.filepath <- file.path(tableDir, paste0("WestfallYoungAdjPvalues_tier2Type3and4.csv"))
  p.df <- read_csv(pval.filepath, show_col_types=FALSE)
  p.df2 <- filter(p.df, mark%in%markName)
  
  p.df <- rbind(p.df1, p.df2)
  p <- as.numeric(p.df[p.df$mark==markName[i], c("p.unadj", "p.FWER", "p.FDR")])
  if(is.na(p[1])){
    unajustedP = sfit$pWald.HRconstant.2sided
    fmt.p <- ifelse(unajustedP<0.001, "< 0.001", paste0("= ", format(unajustedP, digits=2, nsmall=2)))
    subtitle <- paste0("Two-Sided Unadjusted Sieve P ", fmt.p)
  }else{
    fmt.p <-  sapply(p, function(x){ ifelse(x<0.001, "< 0.001", paste0("= ", format(x, digits=2, nsmall=2))) })
    subtitle <- paste0("Two-Sided Unadjusted Sieve P ", fmt.p[1], "\nFWER P ", fmt.p[2], ", Q ", fmt.p[3])
  }
  ylimL <- min(c(sfit$te$TE, -0.4))
  p <- ggplot(sfit, 
              mark=dat1$mark, 
              tx=dat1$tx, 
              xlim=xLim,
              ylim=c(ylimL, 1),
              xtickAt=NULL,
              xtickLab=NULL,
              ytickAt=seq(-1, 1, by=0.2),
              ytickLab=seq(-100, 100, by=20),
              xlab=xLab[i],
              ylab="Vaccine Efficacy (%)",
              axisLabSize=15.5,
              legendLabSize=12,
              txLab=c("Placebo", "Vaccine"),
              jitterFactor=0.1,
              title=title,
              subtitle=subtitle,
              subtitleSize=14,
              estLineSize=2,
              ciLineSize=1.6,
              pointSize=2.1,
              plotHeights=plotHeights)
  print(p)
  dev.off()
}

pVals <- tibble(markName, pHRunity, pHRconstancy)
save(pVals, file=file.path(tabDir, paste0("705_sievePH_pValues_otherQuantMarks.RData")))

