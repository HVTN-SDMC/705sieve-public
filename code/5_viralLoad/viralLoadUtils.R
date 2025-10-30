scientific_10 <- function(x) {
  parse(text=gsub("1e\\+*", " 10^", scales::scientific_format()(x)))
}

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


lm.permuTest <- function(seed = 101, nsim = 2000, formula, data){
  set.seed(seed) ## for reproducibility
  nsim <- nsim
  lmfit <- summary(lm(formula = formula, data = data))
  beta <- lmfit$coefficients["tx", "Estimate"]
  #estimated GM
  allvars <- all.vars(formula)
  coef <- lmfit$coefficients[, "Estimate"]
  W1 <- data.frame("Intercept" = 1,"tx" = 1, data[, allvars[-c(1,2)]])
  W0 <- data.frame("Intercept" = 1,"tx" = 0, data[, allvars[-c(1,2)]])
  
  meanVacs <- mean(as.matrix(W1)%*%coef)
  meanPlacs<- mean(as.matrix(W0)%*%coef)
 
  #permutation tests
  betas_p <- numeric(nsim) ## set aside space for results
  for (i in 1:nsim) {
    pdata <- data
    permTx <- sample(pdata$tx, replace = FALSE)
    pdata$tx <- permTx
    lmfit <- summary(lm(formula = formula, data = pdata))
    betas_p[i] <- lmfit$coefficients["tx", "Estimate"]
  }
  
  #bootstrap confidence interval
  betas_b <- numeric(nsim) ## set aside space for results
  for (i in 1:nsim) {
    permTx <- sample(seq(1, length(data[,1]), 1), replace = TRUE)
    pdata <- data[permTx,]
    lmfit <- summary(lm(formula = formula, data = pdata))
    betas_b[i] <- lmfit$coefficients["tx", "Estimate"]
  }
  
  return(list("P" = sum(abs(betas_p) > abs(beta))/nsim, "CI" = quantile(betas_b, probs = c(0.025, 0.975)), "beta" = beta,
         "meanVac" = meanVacs, "meanPlac" = meanPlacs))
}