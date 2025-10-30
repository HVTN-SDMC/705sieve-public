#Purpose: compare viral load between vaccine and placebo
rm(list = ls(all = TRUE))
library(here)
here::i_am("README.md")
repoDir <- here::here()
dataDir <- file.path(repoDir, "data")
codeDir <- file.path(repoDir, "code/5_viralLoad")
figureDir <- file.path(repoDir, "figures")
tableDir <- file.path(repoDir, "tables")



library(scales)
library(ggpubr)
library(drtmle)
library(SuperLearner)
library(randomForest)
library(np)
library(tidyverse)
library(AFR)
library(MASS)

# load supplement
source(file.path(repoDir, "code/common.R"))
source(file.path(codeDir, "viralLoadUtils.R"))


dat0 <- read.csv(file.path(dataDir, datFile)) %>%
  filter(cohort=="Per-Protocol" & hiv1event == 1)%>%
  dplyr::select(all_of(c("protocol", "subjid", "armdesc", "logrna", "seq.depth.aa", "seq.depth.nt", 
                         "age", "bmi", "riskscoresl", "ind_sa"))) %>%
  mutate(tx = as.numeric(armdesc == "Vaccine"), "vl" = logrna)

dat0$armdesc <- factor(dat0$armdesc, levels = c("Vaccine", "Placebo"))

#vl linear regression
dat <- filter(dat0, !is.na(vl))
W = dplyr::select(dat, all_of(c("age", "bmi","riskscoresl", "ind_sa")))
W$age <- (W$age - mean(W$age, na.rm = TRUE))/sd(W$age)
W$bmi <- (W$bmi - mean(W$bmi, na.rm = TRUE))/sd(W$bmi)
W$riskscoresl <- (W$riskscoresl - mean(W$riskscoresl, na.rm = TRUE))/sd(W$riskscoresl)

lmfit <- lm(dat$vl ~ dat$tx + W$age + W$bmi + W$riskscoresl + W$ind_sa)
lmsummary <- summary(lmfit)
p.vl.firstRNAp <- lmsummary$coefficients["dat$tx", "Pr(>|t|)"]


#amino acid sequencing depth regression
dat <- filter(dat0, !is.na(seq.depth.aa))
W = dplyr::select(dat, all_of(c("age", "bmi","riskscoresl", "ind_sa")))
W$age <- (W$age - mean(W$age, na.rm = TRUE))/sd(W$age)
W$bmi <- (W$bmi - mean(W$bmi, na.rm = TRUE))/sd(W$bmi)
W$riskscoresl <- (W$riskscoresl - mean(W$riskscoresl, na.rm = TRUE))/sd(W$riskscoresl)
dat.pm <- data.frame("seq.depth.aa" = dat$seq.depth.aa, "tx" = dat$tx, "age" = W$age, "bmi" = W$bmi,
                     "riskscoresl" = W$riskscoresl, "ind_sa" = W$ind_sa)
fit <- lm.permuTest(formula = as.formula("seq.depth.aa ~ tx + age + bmi + riskscoresl + ind_sa"),
                                  data = dat.pm)
p.seq.depth.aa.pm <- fit$P


#nt sequencing depth regression
dat <- filter(dat0, !is.na(seq.depth.nt))
W = dplyr::select(dat, all_of(c("age", "bmi","riskscoresl", "ind_sa")))
W$age <- (W$age - mean(W$age, na.rm = TRUE))/sd(W$age)
W$bmi <- (W$bmi - mean(W$bmi, na.rm = TRUE))/sd(W$bmi)
W$riskscoresl <- (W$riskscoresl - mean(W$riskscoresl, na.rm = TRUE))/sd(W$riskscoresl)
p.seq.depth.nt <- lmsummary$coefficients["dat$tx", "Pr(>|t|)"]

dat.pm <- data.frame("seq.depth.nt" = dat$seq.depth.nt, "tx" = dat$tx, "age" = W$age, "bmi" = W$bmi,
                     "riskscoresl" = W$riskscoresl, "ind_sa" = W$ind_sa)
p.seq.depth.nt.pm <- lm.permuTest(formula = as.formula("seq.depth.nt ~ tx + age + bmi + riskscoresl + ind_sa"),
                                  data = dat.pm)$P


my_comparison <- list(c("Vaccine", "Placebo"))
n_vc <- sum(!is.na(dat0$seq.depth.aa[dat$armdesc == "Vaccine"]))
n_pl <- sum(!is.na(dat0$seq.depth.aa[dat$armdesc == "Placebo"]))
set.seed(1)
p <- ggplot() +
  geom_violin(aes(x = armdesc, y = seq.depth.aa, color = armdesc), data = dat) + 
  geom_boxplot(aes(x = armdesc, y = seq.depth.aa, color = armdesc), width = 0.2, data = dat)+
  geom_jitter(aes(x = armdesc, y = seq.depth.aa, color = armdesc),width = 0.2, height=0, data = dat)+
  theme_bw()+
  ylab("Env Amino Acid Sequencing Depth") + 
  xlab("")+
  scale_y_continuous(breaks = seq(0, 600, 100), limits = c(-30, 650), labels = seq(0, 600, 100))+
  geom_segment(aes(x = 2, xend = 2, y = 610, yend = 630))+
  geom_segment(aes(x = 1, xend = 1, y = 610, yend = 630))+
  geom_segment(aes(x = 1, xend = 2, y = 630, yend = 630))+
  annotate (geom = "text", x = 1.5, y = 650, label = paste0("P",format.p(p.seq.depth.aa.pm)), size = 6)+
  annotate (geom = "text", x = 1, y = -30, label = paste0("n = ", n_vc), size = 6)+
  annotate (geom = "text", x = 2, y = -30, label = paste0("n = ", n_pl), size = 6)+
  scale_color_manual(breaks = c("Vaccine", "Placebo"), values = c("darkorange", "darkblue")) +
  guides(color = "none") + 
  theme(legend.title = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        legend.position = c(0.8, 0.95),
        legend.direction = "horizontal",
        legend.text = element_text (size = 16),
        title = element_text(size = 14),
        axis.title.x =  element_text (size = 20),
        axis.title.y =  element_text (size = 20, vjust = +2),
        axis.text.x =  element_text (size = 20),
        axis.text.y =  element_text (size = 20))

ggsave(file.path(figureDir, "AASeqdepthbyTreatment.pdf"), 
       plot = p, height = 7, width = 7)

n_vc <- sum(!is.na(dat0$seq.depth.nt[dat$armdesc == "Vaccine"]))
n_pl <- sum(!is.na(dat0$seq.depth.nt[dat$armdesc == "Placebo"]))
set.seed(1)
p <- ggplot() +
  geom_violin(aes(x = armdesc, y = seq.depth.nt, color = armdesc), data = dat) + 
  geom_boxplot(aes(x = armdesc, y = seq.depth.nt, color = armdesc), width = 0.2, data = dat)+
  geom_jitter(aes(x = armdesc, y = seq.depth.nt, color = armdesc),width = 0.2, height=0, data = dat)+
  theme_bw()+
  ylab("Env Nucleotide Sequencing Depth") + 
  xlab("")+
  scale_y_continuous(breaks = seq(-30, 1000, 200), limits = c(-30, 900), labels = seq(0, 1000, 200))+
  geom_segment(aes(x = 2, xend = 2, y = 850, yend = 870))+
  geom_segment(aes(x = 1, xend = 1, y = 850, yend = 870))+
  geom_segment(aes(x = 1, xend = 2, y = 870, yend = 870))+
  annotate (geom = "text", x = 1.5, y = 900, label = paste0("P", 
                                                            format.p(p.seq.depth.nt.pm)), size = 6)+
  annotate (geom = "text", x = 1, y = -30, label = paste0("n = ", n_vc), size = 6)+
  annotate (geom = "text", x = 2, y = -30, label = paste0("n = ", n_pl), size = 6)+
  scale_color_manual(breaks = c("Vaccine", "Placebo"), values = c("darkorange", "darkblue")) +
  guides(color = "none") + 
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.95),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        legend.direction = "horizontal",
        legend.text = element_text (size = 16),
        title = element_text(size = 14),
        axis.title.x =  element_text (size = 20, vjust = -2),
        axis.title.y =  element_text (size = 20, vjust = +1),
        axis.text.x =  element_text (size = 20),
        axis.text.y =  element_text (size = 20))


ggsave(file.path(figureDir, "NTSeqdepthbyTreatment.pdf"), 
       plot = p, height = 7, width = 7)



###########################################################################################
library(plyr)
viralLoad <- read.csv(file.path(dataDir, "viral_load.csv")) 
viralLoad_preART <- filter(viralLoad, days_since_art >=0 |(is.na(days_since_art) & 
                                                             firstart == "Pre-ART Initiation"))
viralLoad_preART$time <- viralLoad_preART$days_since_dx
viralLoad_preART$subjid <- viralLoad_preART$SUBJID
dat0 <- read.csv(file.path(dataDir, datFile)) %>%
  filter(cohort=="Per-Protocol" & hiv1event == 1)%>%
  dplyr::select(all_of(c("protocol", "subjid", "armdesc", "hiv1fposday"))) %>%
  mutate(tx = as.numeric(armdesc == "Vaccine"))

dat0$armdesc <- factor(dat0$armdesc, levels = c("Vaccine", "Placebo"))
dat0 <- left_join(dat0, viralLoad_preART, by = "subjid")




  
p <- ggplot(data = dat0) +
  geom_point(aes(x = time, y = logrna, group = SUBJID, shape = armdesc, color = armdesc), 
             alpha = 0.8, size = 3, stroke = 1) + 
  geom_line(aes(x = time, y = logrna, group = SUBJID, color = armdesc), 
            alpha = 0.4)+
  xlab("Days since First RNA-Positive Sample")+
  ylab("Viral Load (Copies/ml)")+
  scale_y_continuous(breaks = seq(2, 7, 1), limits = c(1.4, 7.6), labels = scientific_10(10^seq(2, 7, 1)))+
  theme_bw()+
  scale_color_manual(breaks = c("Vaccine", "Placebo"), values = c("darkorange", "darkblue")) + 
  scale_shape_manual(breaks = c("Vaccine", "Placebo"), values = c(5, 21)) + 
  
  theme(legend.title = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        legend.position = c(0.7, 0.95),
        legend.direction = "horizontal",
        legend.text = element_text (size = 16),
        title = element_text(size = 14),
        axis.title.x =  element_text (size = 20),
        axis.title.y =  element_text (size = 20, vjust = +2),
        axis.text.x =  element_text (size = 20),
        axis.text.y =  element_text (size = 20),
        axis.line.x.top = element_blank(), # Remove top x-axis line
        axis.line.y.right = element_blank()) # Remove right y-axis line
ggsave(file.path(figureDir, paste0("longViralLoadbyTreatment.pdf")), 
       plot = p, height = 7, width = 7)
###########################################################################################


#average log10 VL bewteen 2 and 6 weeks
aveViralLoad <- ldply(unique(viralLoad$SUBJID), function(x){
  df <- filter(viralLoad, SUBJID == x & (days_since_art >=0 |(is.na(days_since_art) & 
        firstart == "Pre-ART Initiation") ))
  df <- filter(df, days_since_dx>=14 & days_since_dx <= 42)
  if(length(df$rnacopiesn)>=1){
    vl <- mean(log10(df$rnacopiesn))
    return(c(x, vl))
  }else{
    return(c(x, NA))
  }}
)
colnames(aveViralLoad) <- c("subjid", "ave.vl")
aveViralLoad$subjid <- as.numeric(aveViralLoad$subjid)



dat0 <- read.csv(file.path(dataDir, datFile)) %>%
  filter(cohort=="Per-Protocol" & hiv1event == 1)%>%
  dplyr::select(all_of(c("protocol", "subjid", "armdesc", "logrna", "seq.depth.aa", 
                         "seq.depth.nt", "age", "bmi", "riskscoresl", "ind_sa"))) %>%
  mutate(tx = as.numeric(armdesc == "Vaccine"), "vl" = logrna)

dat0$armdesc <- factor(dat0$armdesc, levels = c("Vaccine", "Placebo"))
dat0 <- left_join(dat0, aveViralLoad , by = "subjid")


dat0$vl.tmp <- as.numeric(dat0[, "ave.vl"])
dat <- filter(dat0, !is.na(vl.tmp) & !is.na(vl))
W = dplyr::select(dat, all_of(c("age", "bmi","riskscoresl", "ind_sa")))
W$age <- (W$age - mean(W$age, na.rm = TRUE))/sd(W$age)
W$bmi <- (W$bmi - mean(W$bmi, na.rm = TRUE))/sd(W$bmi)
W$riskscoresl <- (W$riskscoresl - mean(W$riskscoresl, na.rm = TRUE))/sd(W$riskscoresl)

lmfit <- lm(dat$vl.tmp ~ dat$tx + W$age + W$bmi + W$riskscoresl + W$ind_sa)
lmsummary <- summary(lmfit)
p.vl.ave <- lmsummary$coefficients["dat$tx", "Pr(>|t|)"]

my_comparison <- list(c("Vaccine", "Placebo"))
n_vc <- sum(!is.na(dat$vl.tmp [dat$armdesc == "Vaccine"]))
n_pl <- sum(!is.na(dat$vl.tmp [dat$armdesc == "Placebo"]))



#combining viral load for the first RNA-postive sample and within 2-6 weeks
dat1 <- filter(dat0, !is.na(vl) & !is.na(vl))
dat2 <- filter(dat0, !is.na(vl.tmp) & !is.na(vl))
dat1$armdesc2 <- ifelse(dat1$armdesc == "Vaccine", "1", "2")
dat2$armdesc2 <- ifelse(dat2$armdesc == "Vaccine", "3", "4")

n_vc1 <- sum(!is.na(dat1$vl[dat$armdesc == "Vaccine"]))
n_pl1 <- sum(!is.na(dat1$vl[dat$armdesc == "Placebo"]))
n_vc2 <- sum(!is.na(dat2$vl.tmp [dat$armdesc == "Vaccine"]))
n_pl2 <- sum(!is.na(dat2$vl.tmp [dat$armdesc == "Placebo"]))
set.seed(1)

p <- ggplot() +
  geom_violin(aes(x = armdesc2, y = vl, color = armdesc), data = dat1) + 
  geom_boxplot(aes(x = armdesc2, y = vl, color = armdesc), width = 0.3, data = dat1)+
  geom_jitter(aes(x = armdesc2, y = vl, color = armdesc),width = 0.2, height=0, data = dat1)+
  
  geom_violin(aes(x = armdesc2, y = vl.tmp, color = armdesc), data = dat2) + 
  geom_boxplot(aes(x = armdesc2, y = vl.tmp, color = armdesc), width = 0.3, data = dat2)+
  geom_jitter(aes(x = armdesc2, y = vl.tmp, color = armdesc),width = 0.2, data = dat2)+
  theme_bw()+
  ylab("Viral Load (Copies/ml)") + 
  xlab("")+
  scale_y_continuous(breaks = seq(2, 7, 1), limits = c(1.4, 7.6), labels = scientific_10(10^seq(2, 7, 1)))+
  geom_segment(aes(x = 4, xend = 4, y = 7.127535, yend = 7.318838))+
  geom_segment(aes(x = 3, xend = 3, y = 7.127535, yend = 7.318838))+
  geom_segment(aes(x = 3, xend = 4, y = 7.318838, yend = 7.318838))+
  
  geom_segment(aes(x = 2, xend = 2, y = 7.127535, yend = 7.318838))+
  geom_segment(aes(x = 1, xend = 1, y = 7.127535, yend = 7.318838))+
  geom_segment(aes(x = 1, xend = 2, y = 7.318838, yend = 7.318838))+
  
  annotate (geom = "text", x = 1.5, y = 7.5, label = paste0("First RNA+, P", format.p(p.vl.firstRNAp)), size = 6)+
  annotate (geom = "text", x = 1, y = 1.4, label = paste0("n = ", n_vc1), size = 6)+
  annotate (geom = "text", x = 2, y = 1.4, label = paste0("n = ", n_pl1), size = 6)+
  annotate (geom = "text", x = 3.5, y = 7.5, label = paste0("2-6 Wks, P", format.p(p.vl.ave)), size = 6)+
  annotate (geom = "text", x = 3, y = 1.4, label = paste0("n = ", n_vc2), size = 6)+
  annotate (geom = "text", x = 4, y = 1.4, label = paste0("n = ", n_pl2), size = 6)+
  scale_x_discrete(breaks = c("1", "2", "3", "4"), labels= c("Vaccine", "Placebo", "Vaccine", "Placebo"))+
  scale_color_manual(breaks = c("Vaccine", "Placebo", "Vaccine", "Placebo" ), 
                     values = c("darkorange", "darkblue", "darkorange", "darkblue")) + 
  guides(color = "none")+
  theme(legend.title = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        legend.position = c(0.8, 0.95),
        legend.direction = "horizontal",
        legend.text = element_text (size = 16),
        title = element_text(size = 14),
        axis.title.x =  element_text (size = 20),
        axis.title.y =  element_text (size = 20, vjust = +2),
        axis.text.x =  element_text (size = 20),
        axis.text.y =  element_text (size = 20))

ggsave(file.path(figureDir, paste0("vl_ave_vlbyTreatment.pdf")), 
       plot = p, height = 7, width = 7)





