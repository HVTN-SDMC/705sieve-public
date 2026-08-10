# Purpose: Volcano plot showing the difference in the log HR for the two values of the mark (as the effect size) 
# on the x-axis and –log10 differential VE unadjusted 2-sided p-value on the y-axis. 
# A dashed horizontal line could run through 0.05.
# Author:  Li Li
# Date:    May 9, 2023

# refresh the workspace
rm(list=ls(all=TRUE))

# Setting directory paths -------------------------------------------------
here::i_am("README.md")
repoDir <- here::here()
dataDir <- file.path(repoDir, "data")
codeDir <- file.path(repoDir, "code/2_sieveBinary")
figureDir <- file.path(repoDir, "figures")
tableDir <- file.path(repoDir, "tables")

library(tidyverse)
library(viridis)
source(file.path(repoDir, "code/common.R"))
source(file.path(codeDir,"sieveBinaryUtils.R"))
  
sieveData <- read.csv(file.path(dataDir, dat_file)) %>%
  filter(cohort=="Per-Protocol") %>%
  mutate(tx=as.numeric(armdesc=="Vaccine"), 
         lineageLabel = as.numeric(transmitted.founder.status=="multiple"),
         eventTime = hiv1fposday, 
         eventInd = hiv1event,
         numSequonV5 = as.numeric(num.sequons.v5 >= 2))


posIsAA <- colnames(sieveData)[grepl(".is.",colnames(sieveData))&grepl("hxb2",colnames(sieveData)) & (!grepl("is.sequon.tier1",colnames(sieveData)))]

load(file.path(outputDir, paste0("binaryMarkVE", fileTag,".RData")))


#tier 2 volcano plot
tier2VE <- binaryResultList$tier2
WestfallYoungAdjPvalues <- read.csv(file=file.path(tableDir,"westfallYoung", paste0("WestfallYoungAdjPvalues_tier2Type1and2",fileTag,".csv")))
colnames(WestfallYoungAdjPvalues ) <- c("X", "p.unadj","p.FWER","p.FDR")

tier2Features <- list("posIs" = posIsAA)
for(feature in c("posIs")){
  
  tier2VE.feature <- tier2VE[names(tier2VE) %in% tier2Features[[feature]]]
  
  plot.df <- tibble("mark" = character(),"haplotype"= character() ,"difflogHR" = numeric(), 
                    "neglog10unadjustedp" = numeric(), "unformatedfdr" = numeric(), "unadjustedp" = character(), 
                    "fwer" = character(), "fdr" = character(),"labeltest" = character())
  for(i in 1:length( tier2VE.feature)){
    resulti = tier2VE.feature[[i]]
    VEi = resulti$VEtable$VE
    marki = names(tier2VE.feature)[i]
    markPosi = strsplit(marki, split = "\\.")[[1]][2]
    if(feature == "posIs"){
      haplotypei = strsplit(marki, split = "\\.")[[1]][4]
      haplotypei = paste0(haplotypei, " vs. Not ", haplotypei)
    }else{
      haplotypei = "Match vs. Mismatch"
    }
    difflogHRi = log(1-VEi[1]) -  log(1-VEi[2])
    pvaluesi = filter(WestfallYoungAdjPvalues, X == marki)
    neglog10unadjustedpi = -log10(pvaluesi$p.unadj)
    plot.df <- add_row(.data = plot.df, "mark" = markPosi,"haplotype"= haplotypei ,"difflogHR" = difflogHRi, 
                       "neglog10unadjustedp" = neglog10unadjustedpi, 
                       "unformatedfdr" = pvaluesi$p.FDR,
                       "unadjustedp" = format.p(pvaluesi$p.unadj), "fwer" = format.p0(pvaluesi$p.FWER), 
                       "fdr" = format.p0(pvaluesi$p.FDR),
                       "labeltest" = paste0("Env ", markPosi,"; ", haplotypei,"\n",
                                            "FWER P", format.p2(pvaluesi$p.FWER),"; ",
                                            "Q", format.p2(pvaluesi$p.FDR)))
    
  }
  
  plotLabel <- case_when(feature == "c97za" ~ "Tier 2 Residue Match vs. Mismatch\nClade C C97ZA Boost Env Insert",
                         feature == "posIs" ~ "Tier 2 Residue Scanning")
  p <- ggplot()+
    geom_point(aes(x = difflogHR, y = neglog10unadjustedp, color = round(unformatedfdr,2)), size = 3.5, alpha = 1,
               shape = 21, stroke = 2, fill = "white",
               data = plot.df)+
    # geom_point(aes(x = difflogHR, y = neglog10unadjustedp), size = 4, alpha = 1, color = "red",
    #            data = filter(plot.df, 10^(-neglog10unadjustedp) <= 0.05))+
    #geom_text(aes(x = difflogHR, y = neglog10unadjustedp, label = labeltest), nudge_y = 0.25, 
    #        data = filter(plot.df, 10^(-neglog10unadjustedp) <= 0.05))+
    ggrepel::geom_text_repel(aes(x = difflogHR, y = neglog10unadjustedp, label = labeltest),
                             data = filter(plot.df, 10^(-neglog10unadjustedp) <= 0.05 & unformatedfdr <= 0.2), size = 4.5, 
                             force = 10, seed = 0, lineheight = 0.8, box.padding = 0.5)+
    annotate("segment", x = -0.2, xend = -3, y = -log10(0.0015), yend = -log10(0.0015), linewidth = 1, arrow = arrow(length = unit(0.3, "cm")))+
    annotate("segment", x = 0.2, xend = 3, y = -log10(0.0015), yend = -log10(0.0015), linewidth = 1, arrow = arrow(length = unit(0.3, "cm")))+
    annotate("text", x = -1.7, y = -log10(0.0011), label = "Greater protection against X", size = 5)+
    annotate("text", x = 1.6, y = -log10(0.0011), label = "Greater protection against not X", size = 5)+
    ylab("Differential VE Unadjusted Two-Sided P-value")+
    xlab("Difference in Log Hazard Ratio (Vaccine/Placebo)\nof Residue X vs. Not X")+
    scale_x_continuous(limits = c(-3, 3), breaks = c(-2,-1,0,1, 2), minor_breaks = NULL)+
    scale_y_continuous(limits = c(0, 3), breaks = -log10(c(1, 0.5, 0.2, 0.1,  0.05, 0.01, 0.02, 0.001)),
                       labels = c(1, 0.5, 0.2, 0.1,  0.05, 0.01, 0.02, 0.001), minor_breaks = NULL)+
    #scale_color_steps(name = "Q-value",n.breaks = 5, show.limits = TRUE)+
    #scale_color_viridis(name = "Q-value",n.breaks = 5)+
    #scale_color_steps(name = "Q-value", breaks = c(0.2, 0.4, 0.8), limits = c(0,1), show.limits = TRUE, high = "yellow", low = "red")+
    # scale_color_steps(name = "Q-value",breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    #                   labels = c("0", "0.2", "0.4", "0.6","0.8", "1"),
    #                   limits = c(0,1), high = "green", low = "blue")+
    scale_color_viridis_c(name = "Q-value",breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                          labels = c("0", "0.2", "0.4", "0.6","0.8", "1"),
                          limits = c(0,1), option = "H")+
    theme_bw()+
    theme(plot.margin = unit(c(1,1,1,1), "cm"),
          plot.title = element_text(hjust = 0.5, vjust = 2, size = 18),
          axis.title.x = element_text(size = 17, vjust = -2),
          axis.title.y = element_text(size = 17, vjust = 2),
          axis.text.x = element_text(size = 17, colour = "black"),
          axis.text.y = element_text(size = 17, colour = "black"),
          legend.text=element_text(size=17),
          legend.title = element_text(size=17, vjust =3 ),
          legend.key.size = unit(1.5, "line"))
  
  ggplot2::ggsave(filename = paste0("tier2_", feature, "_Volcano",fileTag, ".pdf"), 
                  plot = p, 
                  path = figureDir,
                  width=0.95*9, height=0.9*7.3)
}






#tier 1 volcano plot
tier1VE <- binaryResultList$tier1
WestfallYoungAdjPvalues <- read.csv(file=file.path(tableDir,"westfallYoung", paste0("WestfallYoungAdjPvalues_tier1Type1to4",fileTag,".csv")))
colnames(WestfallYoungAdjPvalues ) <- c("X", "p.unadj","p.FWER","p.FDR")

tier1Features <- list("posIs" = posIsAA)
for(feature in c("posIs")){
  
  tier1VE.feature <- tier1VE[names(tier1VE) %in% tier1Features[[feature]]]
  plot.df <- tibble("mark" = character(),"haplotype"= character() ,"difflogHR" = numeric(), 
                    "neglog10unadjustedp" = numeric(), "unformatedfdr" = numeric(), "unadjustedp" = character(), 
                    "fwer" = character(), "fdr" = character(),"labeltest" = character())
  for(i in 1:length( tier1VE.feature)){
    resulti = tier1VE.feature[[i]]
    VEi = resulti$VEtable$VE
    marki = names(tier1VE.feature)[i]
    markPosi = strsplit(marki, split = "\\.")[[1]][2]
    if(feature == "posIs"){
      haplotypei = strsplit(marki, split = "\\.")[[1]][4]
      haplotypei = paste0(haplotypei, " vs. Not ", haplotypei)
    }else{
      haplotypei = "Match vs. Mismatch"
    }
    difflogHRi = log(1-VEi[1]) -  log(1-VEi[2])
    pvaluesi = filter(WestfallYoungAdjPvalues, X == marki)
    neglog10unadjustedpi = -log10(pvaluesi$p.unadj)
    plot.df <- add_row(.data = plot.df, "mark" = markPosi,"haplotype"= haplotypei ,"difflogHR" = difflogHRi, 
                       "neglog10unadjustedp" = neglog10unadjustedpi, 
                       "unformatedfdr" = pvaluesi$p.FDR,
                       "unadjustedp" = format.p(pvaluesi$p.unadj), "fwer" = format.p0(pvaluesi$p.FWER), 
                       "fdr" = format.p0(pvaluesi$p.FDR),
                       "labeltest" = paste0("Env ", markPosi,"; ", haplotypei,"\n",
                                            "FWER P", format.p2(pvaluesi$p.FWER),"; ",
                                            "Q", format.p2(pvaluesi$p.FDR)))
    
  }
  
  plotLabel <- case_when(feature == "c97za" ~ "Tier 1 Residue Match vs. Mismatch \nClade C C97ZA Boost Env Insert",
                         feature == "posIs" ~ "Tier 1 Residue Scanning")
  p <- ggplot()+
    geom_point(aes(x = difflogHR, y = neglog10unadjustedp, color = round(unformatedfdr,2)), 
               size = 3.5, alpha = 1, shape = 21, stroke = 2, fill = "white",
               data = plot.df)+
    # geom_point(aes(x = difflogHR, y = neglog10unadjustedp), size = 4, alpha = 1, color = "red",
    #            data = filter(plot.df, 10^(-neglog10unadjustedp) <= 0.05))+
    #geom_text(aes(x = difflogHR, y = neglog10unadjustedp, label = labeltest), nudge_y = 0.25, 
    #        data = filter(plot.df, 10^(-neglog10unadjustedp) <= 0.05))+
    ggrepel::geom_text_repel(aes(x = difflogHR, y = neglog10unadjustedp, label = labeltest),
                             data = filter(plot.df, 10^(-neglog10unadjustedp) <= 0.05 & unformatedfdr< 0.3), size = 5, 
                             force = 5, seed = 0, lineheight = 0.85, box.padding = 0.7,
                             min.segment.length = 1)+
    annotate("segment", x = -0.2, xend = -3, y = -log10(0.0015), yend = -log10(0.0015), linewidth = 1, arrow = arrow(length = unit(0.3, "cm")))+
    annotate("segment", x = 0.2, xend = 3, y = -log10(0.0015), yend = -log10(0.0015), linewidth = 1, arrow = arrow(length = unit(0.3, "cm")))+
    annotate("text", x = -1.7, y = -log10(0.0011), label = "Greater protection against X", size = 5)+
    annotate("text", x = 1.6, y = -log10(0.0011), label = "Greater protection against not X", size = 5)+
    
    ylab("Differential VE Unadjusted Two-Sided P-value")+
    xlab("Difference in Log Hazard Ratio (Vaccine/Placebo)\nof Residue X vs. Not X")+
    
    scale_x_continuous(limits = c(-3, 3), breaks = c(-2,-1,0,1, 2), minor_breaks = NULL)+
    scale_y_continuous(limits = c(0, 3), breaks = -log10(c(1, 0.5, 0.2, 0.1,  0.05, 0.01, 0.02, 0.001)),
                       labels = c(1, 0.5, 0.2, 0.1,  0.05, 0.01, 0.02, 0.001), minor_breaks = NULL)+
    theme_bw()+
    #scale_color_steps(name = "Q-value",breaks = c(0.2, 0.4, 0.6, 0.8), show.limits = TRUE)+
    #scale_color_steps(name = "Q-value", breaks = c(0.2, 0.4, 0.8), limits = c(0,1), show.limits = TRUE, high = "yellow", low = "red")+
    scale_color_viridis_c(name = "Q-value",breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
                          labels = c("0", "0.2", "0.4", "0.6","0.8", "1"),
                          limits = c(0,1), option = "H")+
    theme(plot.margin = unit(c(1,1,1,1), "cm"),
          plot.title = element_text(hjust = 0.5, vjust = 2, size = 18),
          axis.title.x = element_text(size = 17, vjust = -2),
          axis.title.y = element_text(size = 17, vjust = 2),
          axis.text.x = element_text(size = 17, colour = "black"),
          axis.text.y = element_text(size = 17, colour = "black"),
          legend.text=element_text(size=17),
          legend.title = element_text(size=17, vjust =3 ),
          legend.key.size = unit(1.5, "line"))
  
  ggplot2::ggsave(filename = paste0("tier1_", feature, "_Volcano",fileTag, ".pdf"), 
                  plot = p, 
                  path = figureDir,
                  width=0.95*9, height=0.9*7.3)
}