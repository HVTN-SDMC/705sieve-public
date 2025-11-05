#Comparing the structural marks between vaccine and placebo using linear regressions

rm(list = ls(all = TRUE))

here::i_am("README.md")
repoDir <- here::here()
dataDir <- file.path(repoDir, "data")
codeDir <- file.path(repoDir, "code/6_TMscore")
figureDir <- file.path(repoDir, "figures")
tableDir <- file.path(repoDir, "tables")

library(here)
library(tidyverse)
source(file.path(repoDir, "code/common.R"))

source(file.path(codeDir, "utils.R"))
dat <- read.csv(file.path(dataDir, datFile)) %>%
  dplyr::select(all_of(c("subjid", "cohort", "armdesc", "hiv1fposday", "hiv1event", "ind_sa", "age", "bmi", "riskscoresl"))) %>%
  mutate(tx = as.numeric(armdesc == "Vaccine"), eventTime = hiv1fposday, eventInd = hiv1event)

dat <- filter(dat, eventInd == 1 & cohort=="Per-Protocol")

tm <- read.csv(file.path(dataDir, tm_files))
dat1 <- dat %>%
  left_join(tm, by = c("subjid" = "ptid"))

# mark variables to be analyzed
marks <- grep("hxb2", colnames(dat1), value = TRUE)

# keep marks with <=10% NAs among primary endpoints
dat1 <- filter(dat1, eventInd == 1)
marks_keep <- purrr::keep(marks, function(m, df){
  return(mean(is.na(df[, m])) <= 0.1)
}, df = dat1)

marks <- marks_keep
DataSub <- dplyr::select(dat1, all_of(c("eventTime","eventInd","tx","ind_sa","age", "bmi", "riskscoresl", marks)))


df <- sapply(1:length(marks), function(i){
  data1 <- dplyr::select(DataSub,all_of(c("eventTime","eventInd","tx","ind_sa","age", "bmi", "riskscoresl", marks[i])))
  colnames(data1) <- c("eventTime","eventInd","tx","ind_sa", "age", "bmi", "riskscoresl","mark")
  data1 <- filter(data1, eventInd == 1)
  data1$age <- (data1$age - mean(data1$age, na.rm = TRUE))/sd(data1$age)
  data1$bmi <- (data1$bmi - mean(data1$bmi, na.rm = TRUE))/sd(data1$bmi)
  data1$riskscoresl <- (data1$riskscoresl - mean(data1$riskscoresl, na.rm = TRUE))/sd(data1$riskscoresl)
  
  # complete-case analysis, i.e., discard cases with a missing mark
  data1 <- subset(data1, !(eventInd==1 & is.na(mark)))
  data1$mark <- (data1$mark - min(data1$mark))/(max(data1$mark) - min(data1$mark))
  
  lmfit <-lm(mark ~ tx + age + bmi + riskscoresl + ind_sa, data = data1)
  betaest = summary(lmfit)$coefficients["tx", "Estimate"]
  
  return(c(marks[i], betaest))
  
})
df <- data.frame(t(df))
colnames(df) <- c("marks_keep", "beta")
p_df_lb6_v705_A_v1 <- read_csv(file.path(tabDir, "WestfallYoungAdjPvalues_TMscores_lb6_v705_A_v1_lm.csv"))
colnames(p_df_lb6_v705_A_v1)[1] <- "marks_keep"
df <- left_join(df, p_df_lb6_v705_A_v1, by = "marks_keep")
p <- list()
i <- 0

df1 <- df %>%
  mutate(lab = if_else(p.FDR <= 0.2, marks_keep, NA)) %>%
  group_by(beta, p.unadj, p.FDR, p.FWER) %>%
  summarize(lab = paste(unique(lab), collapse = "\n"), .groups = "drop")

for(l in 1:dim(df1)[1]){
  labv <- strsplit(df1$lab[l], split = "\\.")[[1]]
  df1$lab[l] <- paste0(labv[2], "–", labv[4])
}

df1$beta <- as.numeric(df1$beta)
p <- volcano_plot_beta(df1)

ggsave(file.path("figures", paste0("volcano_plot_protomer=A_lb=lb6_lm.pdf")), 
       plot = p, height = 6, width = 8)

