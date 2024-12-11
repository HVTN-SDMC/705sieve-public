# Purpose: Estimation of hazard ratio-based VE by each Tier 1 and Tier 2 Hamming distance 
#          Considers primary endpoints, and the PP cohort is the analysis cohort
#          Hypothesis testing; plotting
#          Time of the first RNA positive sample used as the failure time variable
# Method:  Juraska andn Gilbert (Biometrics, 2013)
#          R package sievePH, version 1.0.3 on Github
# Input:   705 primary endpoints and Hamming distances
# Output:  PDF files, each containing a single plot


# refresh the workspace
rm(list=ls(all=TRUE))

# Setting directory paths -------------------------------------------------
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
sieveData <- read_csv(file.path(datDir, datFile)) %>%
  filter(cohort=="Per-Protocol") %>%
  select(armdesc, eventTime=hiv1fposday, eventInd=hiv1event, ind_sa, starts_with("hdist")) %>%
  mutate(tx=as.numeric(armdesc=="Vaccine"))

######################################################################################################
# Hamming distance variable names
dist <- grep("hdist", colnames(sieveData), value=TRUE)
hdist_ab_no364 <- grep("no364", dist, value = TRUE)
dist <- dist[!dist%in% hdist_ab_no364 ]
# axis labels
xLab <- case_when(dist=="hdist.zspace.mos1.v2" ~ "PC-Weighted Hamming Distance to Mos1 in\nV2 Hotspot in NHP Challenge Model",
                  dist=="hdist.zspace.mos2.v2" ~ "PC-Weighted Hamming Distance to Mos2S in\nV2 Hotspot in NHP Challenge Model",
                  dist=="hdist.zspace.c97za.v2" ~ "PC-Weighted Hamming Distance to C97ZA in\nV2 Hotspot in NHP Challenge Model",
                  dist=="hdist.zspace.mos1.v2_ab" ~ "PC-Weighted Hamming Distance to Mos1 in\nIgG3 V1V2 Correlate Set",
                  dist=="hdist.zspace.mos2.v2_ab" ~ "PC-Weighted Hamming Distance to Mos2S in\nIgG3 V1V2 Correlate Set",
                  dist=="hdist.zspace.c97za.v2_ab" ~ "PC-Weighted Hamming Distance to C97ZA in\nIgG3 V1V2 Correlate Set",
                  dist=="hdist.zspace.mos1.adcc" ~ "PC-Weighted Hamming Distance to Mos1 in\nC1 Epitope of ADCC Correlate in RV144",
                  dist=="hdist.zspace.mos2.adcc" ~ "PC-Weighted Hamming Distance to Mos2S in\nC1 Epitope of ADCC Correlate in RV144",
                  dist=="hdist.zspace.c97za.adcc" ~ "PC-Weighted Hamming Distance to C97ZA in\nC1 Epitope of ADCC Correlate in RV144",
                  dist=="hdist.zspace.mos1.c_ab" ~ "PC-Weighted Hamming Distance to Mos1 in\nClade C bnAb Signatures",
                  dist=="hdist.zspace.mos2.c_ab" ~ "PC-Weighted Hamming Distance to Mos2S in\nClade C bnAb Signatures",
                  dist=="hdist.zspace.c97za.c_ab" ~ "PC-Weighted Hamming Distance to C97ZA in\nClade C bnAb Signatures",
                  dist=="hdist.zspace.mos1.gp41" ~ "gp41 Hamming Distance to Mos1",
                  dist=="hdist.zspace.mos2.gp41" ~ "gp41 Hamming Distance to Mos2S",
                  dist=="hdist.zspace.c97za.gp41" ~ "gp41 Hamming Distance to C97ZA",
                  dist=="hdist.zspace.mos1.gp120" ~ "gp120 Hamming Distance to Mos1",
                  dist=="hdist.zspace.mos2.gp120" ~ "gp120 Hamming Distance to Mos2S",
                  dist=="hdist.zspace.c97za.gp120" ~ "gp120 Hamming Distance to C97ZA",
                  dist=="hdist.zspace.mos1.v5" ~ "V5 Hamming Distance to Mos1",
                  dist=="hdist.zspace.mos2.v5" ~ "V5 Hamming Distance to Mos2S",
                  dist=="hdist.zspace.c97za.v5" ~ "V5 Hamming Distance to C97ZA",
                  dist=="hdist.zspace.1428v1v2" ~ "PC-Weighted Hamming Distance to\ngp70-001428.2.42 V1V2 in IgG3 V1V2 Correlates Set",
                  dist=="hdist.zspace.c97za.hvtn505.cd4bs.antibody" ~ "PC-Weighted Hamming Distance to C97ZA in\nSignature CD4bs Antibody Footprints in HVTN 505",
                  dist=="hdist.zspace.c97za.hvtn505.cd4bs.kmer" ~ "PC-Weighted Hamming Distance to C97ZA in\nCD4bs-Overlapping Signature k-mers",
                  dist=="hdist.zspace.mos1.hvtn505.cd4bs.antibody" ~ "PC-Weighted Hamming Distance to Mos1 in\nSignature CD4bs Antibody Footprints in HVTN 505",
                  dist=="hdist.zspace.mos1.hvtn505.cd4bs.kmer" ~ "PC-Weighted Hamming Distance to Mos1 in\nCD4bs-Overlapping Signature k-mers",
                  dist=="hdist.zspace.mos2.hvtn505.cd4bs.antibody" ~ "PC-Weighted Hamming Distance to Mos2S in\nSignature CD4bs Antibody Footprints in HVTN 505",
                  dist=="hdist.zspace.mos2.hvtn505.cd4bs.kmer" ~ "PC-Weighted Hamming Distance to Mos2S in\nCD4bs-Overlapping Signature k-mers"
                  
)



# strings included in file names
markFileString <- gsub(".", "_", dist, fixed=TRUE)

# initialize vectors of p-values
pHRconstancy <- pHRunity <- NULL

format.p <- function(p, ndigits=2){
  pp <- NULL
  for(i in 1:length(p)){
    if(is.na(p[i])){
      pp[i] <- "--"
    }else{
      if(p[i]<0.001){pp[i] <- " < 0.001"}
      else if (p[i]==1){pp[i] <- "= 1"}
      else{pp[i] <-paste0(" = ",as.character(format(as.numeric(p[i]),digits=ndigits,nsmall=ndigits)))
      if(pp[i]== " = 1.00") {pp[i]= " = 1"}}
    }
  }
  return (pp)
}

# for each Hamming distance
for (i in 1:length(dist)){
  
  dat1 <- select(sieveData, tx, eventTime, eventInd, mark=dist[i], ind_sa)
  
  # convert mark values for non-primary endpoints to NA
  dat1$mark <- ifelse(dat1$eventInd==0, NA, dat1$mark)
  
  # complete-case analysis, i.e., discard cases with a missing mark
  dat1 <- filter(dat1, !(eventInd==1 & is.na(mark)))
  
  # fit the mark-specific HR model
  markRng <- range(dat1$mark, na.rm=TRUE)
  markGrid <- seq(markRng[1], markRng[2], length.out=200)
  fit <- with(dat1, sievePH(eventTime=eventTime, eventInd=eventInd, mark=mark, tx=tx, strata=ind_sa))
  sfit <- summary(fit, markGrid=markGrid, sieveAlternative="oneSided")
  
  # 2-sided Wald test of {H0: VE(v)=0 for all v}
  pHRunity <- c(pHRunity, sfit$pWald.HRunity.2sided)
  
  # Double 1-sided Wald test of {H0: VE(v)=VE for all v} against {H1: VE(v) varies with v}
  
  pHRconstancy <- c(pHRconstancy, min(2*sfit$pWald.HRconstant.1sided, 1))
  
  if(dist[i] %in% c("hdist.zspace.c97za.c_ab","hdist.zspace.mos1.c_ab","hdist.zspace.mos2.c_ab" )){
    xtickAt = c(0, 2, 4, 6,8)
    xtickLab = c(0, 2, 4, 6,8)
    xLim <- c(0,8)
  }else if(dist[i] %in% c("hdist.zspace.c97za.hvtn505.cd4bs.kmer")){
    xtickAt = c(0, 2, 4, 6,8,10)
    xtickLab = c(0, 2, 4, 6,8,10)
    xLim = c(0,10)
  }else if (dist[i] %in% c("hdist.zspace.mos2.hvtn505.cd4bs.kmer")){
    xtickAt = c(2, 4, 6,8,10)
    xtickLab = c(2, 4, 6,8,10)
    xLim = c(2,10)
  }else if (dist[i] %in% c("hdist.zspace.mos1.hvtn505.cd4bs.kmer")){
    xtickAt = c(4, 6,8,10)
    xtickLab = c(4, 6,8,10)
    xLim = c(4,10)
  }else if(dist[i] %in% c("hdist.zspace.c97za.v2", "hdist.zspace.mos1.v2", "hdist.zspace.mos2.v2")){
    xtickAt = c(0, 5, 10, 15)
    xtickLab = c(0, 5, 10, 15)
    xLim = c(0,15)
  }else if(dist[i] %in% c("hdist.zspace.c97za.v2_ab", "hdist.zspace.mos1.v2_ab", "hdist.zspace.mos2.v2_ab")){
    xtickAt = c(0, 5, 10, 15)
    xtickLab = c(0, 5, 10, 15)
    xLim = c(0,15)
    
  }else if(dist[i] %in% c("hdist.zspace.c97za.adcc", "hdist.zspace.mos1.adcc", "hdist.zspace.mos2.adcc")){
    xtickAt = c(0, 1, 2, 3,4)
    xtickLab = c(0, 1, 2, 3,4)
    xLim = c(0,4)
  }else if(dist[i] %in% c("hdist.zspace.c97za.hvtn505.cd4bs.antibody", "hdist.zspace.mos1.hvtn505.cd4bs.antibody", "hdist.zspace.mos2.hvtn505.cd4bs.antibody")){
    xtickAt = c(5, 10, 15,20)
    xtickLab = c(5, 10, 15,20)
    xLim = c(5,20)
  }else{
    xLim <- range(dat1$mark, na.rm=TRUE)
    xtickAt = NULL
    xtickLab = NULL
  }
  rg = range(dat1$mark, na.rm=TRUE)
  xLim = c(min(xLim[1], rg[1]), max(xLim[2], rg[2]))
  
  
  
  pdf(file.path(figDir, paste0("705_sievePH_VE_", markFileString[i], "_PP.pdf")), width=0.9*7, height=0.9*6.3)
  
  plotHeights <- c(0.32, 0.68)
  ggsave.width <- 0.7 * 7.3
  ggsave.height <- 0.7 * 6.5
  title <- NULL
  pval.filepath <- file.path(tableDir, paste0("WestfallYoungAdjPvalues_tier1Type5to7.csv"))
  p.df <- read_csv(pval.filepath, show_col_types=FALSE)
  p.df1 <- filter(p.df, mark%in%dist)
  pval.filepath <- file.path(tableDir, paste0("WestfallYoungAdjPvalues_tier2Type3and4.csv"))
  p.df <- read_csv(pval.filepath, show_col_types=FALSE)
  p.df2 <- filter(p.df, mark%in%dist)
  
  p.df <- rbind(p.df1, p.df2)
  p <- as.numeric(p.df[p.df$mark==dist[i], c("p.unadj", "p.FWER", "p.FDR")])
  if(is.na(p[1])){
    unajustedP = min(2*sfit$pWald.HRconstant.1sided,1)
    fmt.p <- format.p(unajustedP)
    subtitle <- paste0("Double One-Sided Unadjusted Sieve P ", fmt.p)
  }else{
    fmt.p <-  sapply(p, function(x){ format.p(x)})
    subtitle <- paste0("Double One-Sided Unadjusted Sieve P ", fmt.p[1], "\nFWER P ", fmt.p[2], ", Q ", fmt.p[3])
  }
  ylimL <- -0.4
  
  p <- ggplot.summary.sievePH(sfit, 
                              mark=dat1$mark, 
                              tx=dat1$tx, 
                              xlim=xLim,
                              ylim=c(ylimL, 1),
                              xtickAt=xtickAt,
                              xtickLab=xtickLab,
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

pVals <- tibble(dist, pHRunity, pHRconstancy)
save(pVals, file=file.path(tabDir, paste0("705_sievePH_pValues_hdist.RData")))


