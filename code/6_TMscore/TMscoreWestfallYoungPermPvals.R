#  Westfall and Young permutation-based multiplicity adjustment for TM-score mark variables

#Comparing the structural marks between vaccine and placebo using linear regressions

rm(list = ls(all = TRUE))

here::i_am("README.md")
repoDir <- here::here()
dataDir <- file.path(repoDir, "data")
codeDir <- file.path(repoDir, "code/6_TMscore")
figureDir <- file.path(repoDir, "figures")
tableDir <- file.path(repoDir, "tables")

library(tidyverse)

source(file.path(repoDir, "code/common.R"))
source(file.path(repoDir, "code/6_TMscore/p.adj.perm2.R"))


##########################################################################################################################
#Comparing the structural marks between vaccine and placebo using linear regressions
dat <- read.csv(file.path(dat_path, dat_file)) %>%
  dplyr::select(all_of(c("cohort","subjid", "armdesc", "hiv1fposday", "hiv1event", "ind_sa", "age", "bmi", "riskscoresl"))) %>%
  mutate(tx = as.numeric(armdesc == "Vaccine"), eventTime = hiv1fposday, eventInd = hiv1event)

dat <- filter(dat, eventInd == 1 & cohort=="Per-Protocol")

tm <- read.csv(file.path(dat_path, tm_files))
dat1 <- dat %>%
  left_join(tm, by = c("subjid" = "ptid"))

# mark variables to be analyzed
marks <- grep("hxb2", colnames(dat1), value = TRUE)

# keep marks with <=10% NAs among primary endpoints
edat1 <- filter(dat1, eventInd == 1)
marks_keep <- purrr::keep(marks, function(m, df){
  return(mean(is.na(df[, m])) <= 0.1)
}, df = edat1)

nPerm <- 1000

lm.permuTest <- function(nsim = 1000, formula, data){
  nsim <- nsim
  n = dim(data)[1]
  X <- cbind(rep(1, n), data$tx, data$age, data$bmi, data$riskscoresl, data$ind_sa)
  Y <- data$mark
  betav <- solve(t(X)%*%X)%*%t(X)%*%Y
  beta <- betav[2]
  betas <- numeric(nsim) ## set aside space for results
  for (i in 1:nsim) {
    pdata <- data
    permTx <- sample(pdata$tx, replace = FALSE)
    X <- cbind(rep(1, n), permTx, data$age, data$bmi, data$riskscoresl, data$ind_sa)
    Y <- data$mark
    betav <- solve(t(X)%*%X)%*%t(X)%*%Y
    betas[i] <- betav[2]
  }
  return(sum(abs(betas) > abs(beta))/nsim)
}

# Compute p-values from data sets with resampled marks --------------------
marks <- marks_keep
DataSub <- dplyr::select(dat1, all_of(c("eventTime","eventInd","tx","ind_sa", "age", "bmi", "riskscoresl", marks)))
# get the p-values for individual permutations of the observed marks
# the first vector in 'pvals' stores unadjusted p-values based on original data

pvals <- lapply(1:(nPerm + 1), function(seed){
  set.seed(seed)
  print(seed)
  # permute the marks observed or missing in cases
  idx <- sample(1:sum(DataSub$eventInd))
  
  # 'pvals1' is a vector of p-values (one for each mark variable) for the single permutation
  pvals1 <- sapply(1:length(marks), function(i){
    data1 <- dplyr::select(DataSub,all_of(c("eventTime","eventInd","tx","ind_sa","age", "bmi", "riskscoresl", marks[i])))
    colnames(data1) <- c("eventTime","eventInd","tx","ind_sa", "age", "bmi", "riskscoresl","mark")
    data1 <- filter(data1, eventInd == 1)
    data1$age <- (data1$age - mean(data1$age, na.rm = TRUE))/sd(data1$age)
    data1$bmi <- (data1$bmi - mean(data1$bmi, na.rm = TRUE))/sd(data1$bmi)
    data1$riskscoresl <- (data1$riskscoresl - mean(data1$riskscoresl, na.rm = TRUE))/sd(data1$riskscoresl)
    
     # apply the permutation
    if (seed>1){
      data1$tx[data1$eventInd==1] <- data1$tx[data1$eventInd==1][idx]  
    }
    
    # complete-case analysis, i.e., discard cases with a missing mark
    data1 <- subset(data1, !(eventInd==1 & is.na(mark)))
    data1$mark <- (data1$mark - min(data1$mark))/(max(data1$mark) - min(data1$mark))
    
    dat.pm <- data.frame("mark" = data1$mark, "tx" = data1$tx, "age" = data1$age, "bmi" = data1$bmi,
                         "riskscoresl" = data1$riskscoresl, "ind_sa" = data1$ind_sa)
    

    p.mark.pm <- lm.permuTest(formula = as.formula("mark ~ tx + age + bmi + riskscoresl + ind_sa"),
                              data = dat.pm)
   
    return(p.mark.pm)
  })
  
  names(pvals1) <- marks
  return(pvals1)
})

pvals <- do.call(rbind, pvals)
# Apply Westfall and Young (1993) to obtain adjusted p-values -------------
pvals.adj <- p.adj.perm2(p.unadj=pvals[1, ], p.perms=pvals[-1, ], alpha=1)
mark <- as.vector(row.names(pvals.adj))
npvals.adj <- data.frame(mark,pvals.adj )
write.csv(npvals.adj, file=file.path(tableDir, paste0("WestfallYoungAdjPvalues_TMscores_lb6_v705_A_v1_lm.csv")), row.names=FALSE)
save(pvals, file=file.path(tableDir, paste0("WestfallYoungPermPvalues_TMscores_lb6_v705_A_v1_lm.RData")))



