library(dplyr)
library(DT)
library(purrr)
library(ggplot2)
library(ggpubr)
library(survival)
library(deconvolveR)
library(REBayes)
library(MASS)


#' Construct a prior on event probabilities.
#'
#' This function builds a prior distribution for the binomial probability
#' \eqn{\theta} based on observed counts \code{K} and denominators \code{M}.
#' Depending on \code{prior_type}, the prior is either:
#' \itemize{
#'   \item \code{"raw"}: empirical point-mass at the observed ratios \code{K / M};
#'   \item \code{"spline"}: a smoothed deconvolution estimate via \code{deconv_fast};
#'   \item \code{"npmle"}: a nonparametric MLE via \code{REBayes::Bmix}.
#' }
#'
#' @param K Integer (or numeric) vector of observed counts (e.g., number of events).
#' @param M Integer (or numeric) vector of denominators (e.g., sequencing depth / trials),
#'   the same length as \code{K}.
#' @param prior_type Character string, one of \code{"raw"}, \code{"spline"} or
#'   \code{"npmle"}, controlling how the prior is estimated.
#' @param degree_spline Integer degree of the spline passed to \code{deconv_fast}
#'   when \code{prior_type = "spline"}.
#' @param reg_param Nonnegative regularization parameter \code{c0} passed to
#'   \code{deconv_fast} when \code{prior_type = "spline"}.
#' @param prior_grid Numeric vector giving the support grid \code{v} used by
#'   \code{REBayes::Bmix} when \code{prior_type = "npmle"}.
#' @return A data frame with two columns:
#' \describe{
#'   \item{theta}{Grid of probability values.}
#'   \item{g}{Prior mass (or density) at each value of \code{theta}.}
#' }
obtain_prior <- function(K, M, prior_type = "spline", degree_spline = 15, 
                         reg_param = 1, prior_grid = seq(0, 1, 0.005)){
  
  if (length(K) != length(M)) {
    stop(sprintf("K and M must have the same length. Got lengths %d and %d.",
                 length(K), length(M)))
  }
  if (any(!is.finite(K)) || any(!is.finite(M))) {
    stop("K and M must be finite numeric vectors (no NA/NaN/Inf).")
  }
  if (any(M <= 0)) {
    stop("All M must be > 0 (cannot divide by zero or use non-positive trials).")
  }
  
  # ---- prior_type validation with friendly message ----
  allowed <- c("raw", "spline", "npmle")
  if (!prior_type %in% allowed) {
    stop(sprintf("Invalid prior_type: '%s'. Must be one of: %s.",
                 prior_type, paste(allowed, collapse = ", ")))
  }
  
  # helper: empirical point-mass density on unique grid points of K/M
  dirac_density <- function(data) {
    unique_points <- unique(data) 
    density_values <- sapply(unique_points, function(pt) {
      sum(data == pt) / length(data)
    })
    data.frame(theta = unique_points, g = density_values)
  }
  
  if (prior_type == "raw"){
    prior <- dirac_density(K / M)
  }
  else if (prior_type == "spline"){
    tau <- seq(0, 1, 0.002)
    df_deconv <- data.frame(n = M, 
                            s = K)
    deconv_res <- deconv_fast(tau = tau, X = df_deconv, family = "Binomial", 
                         c0 = reg_param, pDegree = degree_spline)
    prior <- data.frame(deconv_res$stats)
    
  }
  else if (prior_type == "npmle"){
    npmle_res <- Bmix(K, M, v = prior_grid)
    
    prior <- data.frame(theta = npmle_res$x, 
                        g = npmle_res$y)
  }
  prior
}


#' Posterior probability that the event rate exceeds a threshold.
#'
#' Given binomial counts \code{K} and denominators \code{M}, this function
#' first estimates a prior on the true event probability \eqn{\theta} using
#' \code{\link{obtain_prior}}, then computes, for each observation \code{i},
#' the posterior probability \eqn{P(\theta_i \ge \text{threshold} \mid K_i, M_i)}.
#'
#' @param K Integer (or numeric) vector of observed counts.
#' @param M Integer (or numeric) vector of denominators (e.g., sequencing depth),
#'   the same length as \code{K}.
#' @param threshold Numeric scalar giving the cutpoint for \eqn{\theta}; the
#'   function returns posterior probabilities that \eqn{\theta \ge \text{threshold}}.
#' @param prior_type Character string passed to \code{\link{obtain_prior}},
#'   one of \code{"raw"}, \code{"spline"}, or \code{"npmle"}.
#' @param degree_spline Integer degree for the spline prior when
#'   \code{prior_type = "spline"}.
#' @param reg_param Regularization parameter \code{c0} for \code{deconv_fast}
#'   when \code{prior_type = "spline"}.
#' @param prior_grid Numeric vector specifying the support grid for the
#'   NPMLE prior when \code{prior_type = "npmle"}.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{K}{Observed counts.}
#'   \item{sequence_depth}{Denominators \code{M}.}
#'   \item{peen}{Posterior probability that \eqn{\theta \ge \text{threshold}}.}
#'   \item{pnul}{Complementary probability \code{1 - peen}.}
#' }
estimate_posterior_probs <- function(K, M, threshold = 0.01, 
                                     prior_type = "spline", 
                                     degree_spline = 15, 
                                     reg_param = 1, 
                                     prior_grid = seq(0, 1, 0.005)){
  
  prior <- obtain_prior(K, M, prior_type, degree_spline, reg_param, 
                        prior_grid)
  
  p_gt_threshold <- sapply(1:length(M), function(i) {
    
    likelihood <- dbinom(K[i], size = M[i], prob = prior$theta)
    posterior <- prior$g * likelihood
    
    posterior <- posterior / sum(posterior) 
    
    # Identify grid points where Q >= threshold
    q_greater <- prior$theta >= threshold
    
    # Sum posterior probabilities for Q > threshold
    sum(posterior[q_greater])
  })
  
  data.frame(
    K = K, 
    sequence_depth = M, 
    peen = p_gt_threshold, 
    pnul = 1 - p_gt_threshold
  )
}


#' Add classification probabilities to a data frame by groups.
#'
#' For each distinct value of \code{group} in \code{df}, this function:
#' \enumerate{
#'   \item Subsets to rows with non-missing \code{K} in that group;
#'   \item Calls \code{\link{estimate_posterior_probs}} on \code{K} and
#'         \code{sequence_depth};
#'   \item Merges the resulting posterior probabilities back to the original
#'         data set by \code{K}, \code{sequence_depth}, and \code{group}.
#' }
#' For censored observations (\code{allstatus == 0}), both classification
#' probabilities \code{peen} and \code{pnul} are set to 0.
#'
#' @param df A data frame containing at least the columns:
#'   \code{K} (counts), \code{sequence_depth} (denominators),
#'   \code{group} (grouping factor or id), and \code{allstatus}
#'   (event/censoring indicator).
#' @param thresh Numeric threshold passed to
#'   \code{\link{estimate_posterior_probs}} as \code{threshold}.
#' @param prior_type Character string specifying how the prior is obtained;
#'   passed to \code{\link{estimate_posterior_probs}} (and ultimately
#'   \code{\link{obtain_prior}}).
#' @param degree_spline Either a single integer, used for all groups when
#'   \code{prior_type = "spline"}, or a vector of integers, one per group.
#' @param reg_param Regularization parameter for the spline prior when
#'   \code{prior_type = "spline"}.
#' @param prior_grid Numeric grid of support values for the NPMLE prior when
#'   \code{prior_type = "npmle"}.
#'
#' @return A copy of \code{df} with two additional columns:
#' \describe{
#'   \item{peen}{Posterior probability that \eqn{\theta \ge \text{thresh}}.}
#'   \item{pnul}{Complementary probability \code{1 - peen}.}
#' }
#'
add_classification_probs <- function(df, 
                                     thresh = 0.01, 
                                     prior_type = "spline",
                                     degree_spline = 15, 
                                     reg_param = 1, 
                                     prior_grid = seq(0, 1, 0.005)){
  
  df_to_merge <- data.frame()
  
  idx <- 1
  
  for (group_level in unique(df$group)){
    
    if (length(degree_spline) == 1){
      degree_spline_temp <- degree_spline
    }
    else{
      degree_spline_temp <- degree_spline[idx]
      idx <- idx + 1
    }
    
    only_available_K_df <- df %>% 
      filter(!is.na(K)) %>% 
      filter(group == group_level)
    
    df_to_merge_temp <- estimate_posterior_probs(K = only_available_K_df$K, 
                                                 M = only_available_K_df$sequence_depth, 
                                                 prior_type = prior_type, 
                                                 threshold = thresh, 
                                                 degree_spline = degree_spline_temp,
                                                 reg_param = reg_param,
                                                 prior_grid = prior_grid) %>% 
      distinct(K, sequence_depth, peen, pnul) %>% 
      mutate(group = group_level)
    

    df_to_merge <- rbind(df_to_merge, df_to_merge_temp)
    
  }
  
  df %>% 
    merge(df_to_merge, all.x = T) %>% 
    
    # set both classification probabilities as 0 if censored
    mutate(peen = ifelse(allstatus == 0, 0, peen), 
           pnul = ifelse(allstatus == 0, 0, pnul))
}


#' Stratified partial log-likelihood with classification weights.
#'
#' This function evaluates the (negative) weighted partial log-likelihood
#' for a two-cause proportional hazards model, with cause-specific hazards
#' parameterized by \code{rho} and \code{phi}, and observation-specific
#' auxiliary weights \code{pnul} and \code{peen} giving the probabilities
#' of each cause.
#'
#' The data are assumed to be sorted in increasing event time within each
#' stratum (so that reverse cumulative sums over the sorted rows define the
#' risk sets). The design matrix is constructed from the covariates in
#' \code{X_vars} and \code{Z_vars}, and the parameter vector \code{par}
#' is mapped to the cause-specific linear predictors via the constrained
#' layout used in the original implementation.
#'
#' @param par Numeric vector of free parameters. Its length should be
#'   \code{length(X_vars) + length(Z_vars)}, corresponding to the non-zero
#'   entries in the constrained coefficient vectors for the two causes.
#' @param data A data frame containing:
#'   \itemize{
#'     \item the covariates named in \code{X_vars} and \code{Z_vars};
#'     \item the classification weights \code{pnul} and \code{peen};
#'     \item variables defining strata, if \code{strata_var} is not \code{NULL}.
#'   }
#'   Rows must be ordered by increasing event time within strata before calling.
#' @param X_vars Character vector of column names in \code{data} to be used
#'   as \eqn{X}-covariates (for the first cause).
#' @param Z_vars Character vector of column names in \code{data} to be used
#'   as \eqn{Z}-covariates (for the second cause).
#' @param strata_var Optional character vector of column name(s) in \code{data}
#'   that define the strata. If \code{NULL}, all observations are treated as a
#'   single stratum.
#' @param pre Optional list of precomputed quantities, typically created by
#'   \code{\link{build_pre}}. When supplied, the function uses
#'   \code{pre$datamat}, \code{pre$strata}, \code{pre$nX}, and \code{pre$nZ}
#'   instead of reconstructing them from \code{data}, \code{X_vars},
#'   \code{Z_vars}, and \code{strata_var}.
#'
#' @return A single numeric value: the negative weighted partial log-likelihood.
#'   This is suitable for use as an objective function in optimizers such as
#'   \code{\link[stats]{optim}}.
#'   
loglik_strat <- function(par, data, X_vars, Z_vars, strata_var = NULL, pre = NULL) {

  if (!is.null(pre)) {
    # use precomputed pieces
    nX      <- pre$nX
    nZ      <- pre$nZ
    datamat <- pre$datamat
    strata  <- pre$strata
    n       <- nrow(datamat)
  } else {
    # original behavior
    nX <- length(X_vars)
    nZ <- length(Z_vars)
    datamat <- as.matrix(data[, c(X_vars, Z_vars), drop = FALSE])
    n <- nrow(datamat)
    strata <- if (is.null(strata_var)) factor(rep(1L, n)) else {
      interaction(data[, strata_var, drop = FALSE], drop = TRUE, lex.order = TRUE)
    }
  }
  
  ## Parameter layout (unchanged)
  if (nX == 0) {
    rho <- if (nZ == 0) numeric(0) else rep(0, nZ)
  } else {
    rho <- c(if (nZ == 0) NULL else rep(0, nZ), par[1:nX])
  }
  if (nZ == 0) {
    phi <- if (nX == 0) numeric(0) else rep(0, nX)
  } else {
    phi <- c(par[(nX + 1):(nX + nZ)], if (nX == 0) NULL else rep(0, nX))
  }
  
  ## Linear predictors and hazards (unchanged)
  eta_rho <- if (length(rho)) as.vector(datamat %*% rho) else rep(0, n)
  eta_phi <- if (length(phi)) as.vector(datamat %*% phi) else rep(0, n)
  haz0    <- exp(eta_rho)
  haz1    <- exp(eta_phi)
  
  ## reverse cumsum within strata (unchanged)
  rev_cumsum_by <- function(x, g) {
    ave(x, g, FUN = function(v) rev(cumsum(rev(v))))
  }
  
  ## Risk-set totals (unchanged)
  risk0 <- rev_cumsum_by(haz0, strata)
  risk1 <- rev_cumsum_by(haz1, strata)
  
  ## Auxiliary weights (unchanged)
  pi0  <- data$pnul
  pi1  <- data$peen
  fail <- (pi0 + pi1) > 0
  
  ## Stratified partial log-likelihood (unchanged)
  ll_vec <- pi0 * (log(haz0) - log(risk0)) +
    pi1 * (log(haz1) - log(risk1))
  
  -sum(ll_vec[fail])
}


# helper function 
build_pre <- function(data, X_vars, Z_vars, strata_var = NULL) {
  datamat <- as.matrix(data[, c(X_vars, Z_vars), drop = FALSE])
  n       <- nrow(datamat)
  strata  <- if (is.null(strata_var)) factor(rep(1L, n)) else {
    interaction(data[, strata_var, drop = FALSE], drop = TRUE, lex.order = TRUE)
  }
  list(
    datamat = datamat,
    strata  = strata,
    nX = length(X_vars),
    nZ = length(Z_vars)
  )
}


#' Run a proportional hazards analysis using posterior classification probabilities.
#'
#' This function estimates posterior classification probabilities for each
#' observation based on auxiliary variables, incorporates those probabilities
#' into a weighted partial likelihood, and fits a proportional hazards model
#' using the optimization routine in \code{optim()}.
#'
#' Steps performed:
#' \enumerate{
#'   \item Construct a grouping variable from \code{cl_probs_vars}.
#'   \item Compute posterior classification probabilities within each group via
#'         \code{\link{add_classification_probs}}.
#'   \item Order the data by event time (censorings pushed to the end on ties).
#'   \item Build the design matrices using \code{\link{build_pre}}.
#'   \item Maximize the weighted stratified partial likelihood
#'         \code{\link{loglik_strat}} using \code{optim()}.
#' }
#'
#' @param df Data frame containing event times, event indicators, covariates,
#'   and the variables needed to compute classification probabilities.
#' @param cl_probs_vars Character vector of column names whose combined values
#'   define the auxiliary grouping used to estimate priors.
#' @param covariates Character vector of covariate names used in both
#'   cause-specific linear predictors.
#' @param event_var Name of the event/censoring indicator column.
#' @param time_var Name of the observed time variable.
#' @param prior_type Character string passed through to
#'   \code{\link{add_classification_probs}} → \code{\link{estimate_posterior_probs}}
#'   → \code{\link{obtain_prior}}.
#' @param stratification_var Optional column name(s) defining strata. If
#'   \code{NULL}, all observations belong to one stratum.
#' @param degree_spline Spline degree for the prior when
#'   \code{prior_type = "spline"}.
#' @param reg_param Regularization parameter for the spline prior.
#' @param prior_grid Support grid for the nonparametric MLE prior.
#' @param par_init Initial parameter vector passed to \code{optim()}.
#' @param thresh Threshold used when computing posterior classification
#'   probabilities.
#'
#' @return A numeric vector of estimated parameters (the optimized output from
#'   \code{optim()}).
#'   
run_corrected_analysis <- function(df, cl_probs_vars, covariates, 
                                   event_var = "allstatus",
                                   time_var = "observed_time",
                                   prior_type = "spline",
                                   stratification_var = NULL,
                                   degree_spline = 15, 
                                   reg_param = 1, 
                                   prior_grid = seq(0, 1, 0.005), 
                                   par_init = c(0.69, -0.05),
                                   thresh = 0.01){
  
  # create group variable in df 
  df$group <- as.integer(interaction(df[cl_probs_vars], drop = TRUE))
  
  # add classification probabilities by group 
  df_new <- add_classification_probs(df, 
                                     thresh = thresh, 
                                     prior_type = prior_type,
                                     degree_spline = degree_spline, 
                                     reg_param = reg_param, 
                                     prior_grid = prior_grid)
  
  # order by event/censor  time 
  timeorder <- function(data){
    data$times <- data[,time_var] + 1e-12 * (1 - data[,event_var]) # censorings last
    return(data[order(data$times), ])
  }
  
  df_new <- timeorder(df_new)
  
  pre <- build_pre(df_new, X_vars = covariates, Z_vars = covariates,
                   strata_var = stratification_var)
  
  opt_res <- optim(par = par_init,
                   fn  = loglik_strat,
                   method = "BFGS",
                   data = df_new, X_vars = covariates, Z_vars = covariates,
                   strata_var = stratification_var,
                   pre = pre) 
  
  par_res <- opt_res$par
  
  par_res
}


#' Bootstrap for the classification-weighted analysis.
#'
#' This function repeatedly resamples the rows of \code{df} with replacement,
#' runs \code{\link{run_corrected_analysis}} on each bootstrap sample, and
#' stores the resulting parameter estimates.
#'
#' @param df Data frame used as the original dataset for bootstrapping.
#' @param cl_probs_vars Character vector of column names defining the grouping
#'   used for estimating classification probabilities.
#' @param covariates Character vector of covariate names used in the model.
#' @param event_var Name of the event/censoring indicator column.
#' @param time_var Name of the observed time variable.
#' @param stratification_var Optional column name(s) defining strata.
#' @param boots Integer number of bootstrap resamples.
#' @param prior_type Character string specifying how the prior is obtained;
#'   passed through to \code{\link{run_corrected_analysis}}.
#' @param degree_spline Spline degree for the prior when
#'   \code{prior_type = "spline"}.
#' @param reg_param Regularization parameter for the spline prior.
#' @param prior_grid Support grid for the nonparametric MLE prior.
#' @param par_init Initial parameter vector passed to
#'   \code{\link{run_corrected_analysis}} and ultimately \code{optim()}.
#' @param thresh Threshold used when computing posterior classification
#'   probabilities.
#'
#' @return A numeric matrix with \code{boots} rows and
#'   \code{length(covariates) * 2} columns, where each row contains one set of
#'   bootstrap parameter estimates. Rows corresponding to failed optimizations
#'   (if any) remain \code{NA}.
#'
run_boot_corrected_analysis <- function(df, cl_probs_vars, covariates, 
                                        event_var = "allstatus",
                                        time_var = "observed_time",
                                        stratification_var = NULL, 
                                        boots = 100,
                                        prior_type = "spline",
                                        degree_spline = 15, 
                                        reg_param = 1, 
                                        prior_grid = seq(0, 1, 0.005), 
                                        par_init = c(0.69, -0.05), 
                                        thresh = 0.01){

  boot_estimates <- matrix(NA, nrow = boots, ncol = length(covariates) * 2)
  
  for (b in 1:boots) {
    boot_df <- df[sample(1:nrow(df), replace = TRUE), ]
    
    try({
    res_temp <- run_corrected_analysis(boot_df, 
                                       cl_probs_vars = cl_probs_vars, 
                                       covariates = covariates, 
                                       event_var = event_var,
                                       time_var = time_var,
                                       stratification_var = stratification_var,
                                       prior_type = prior_type,
                                       degree_spline = degree_spline, 
                                       reg_param = reg_param, 
                                       prior_grid = prior_grid, 
                                       par_init = par_init, 
                                       thresh = thresh)
    boot_estimates[b, ] <- res_temp
    
    })
  }
  
  boot_estimates
}


#' Wald-type joint test for the null hypothesis of no vaccine effectiveness.
#'
#' This function performs a two-parameter Wald test of
#' \deqn{H_0: (\beta_0, \beta_1) = (0, 0)}
#' using the bootstrap covariance matrix of the parameter estimates.
#'
#' @param overall_estimates Numeric vector of parameter estimates from the
#'   full dataset (typically output of \code{run_corrected_analysis}).
#' @param boot_estimates Numeric matrix where each row contains bootstrap
#'   parameter estimates (typically from \code{run_boot_corrected_analysis}).
#' @param col_beta0 Integer column index of \code{boot_estimates} containing
#'   the estimates for \eqn{\beta_0}.
#' @param col_beta1 Integer column index of \code{boot_estimates} containing
#'   the estimates for \eqn{\beta_1}.
#'
#' @return A single numeric value: the p-value of the Wald chi-square test
#'   with 2 degrees of freedom.
#'
hypothesis_test_no_VE <- function(overall_estimates, 
                                  boot_estimates, 
                                  col_beta0 = 1, 
                                  col_beta1 = 2){
  
  beta0_est <- overall_estimates[col_beta0]
  beta1_est <- overall_estimates[col_beta1]
  
  beta0_boot <- boot_estimates[,col_beta0]
  beta1_boot <- boot_estimates[,col_beta1]
  
  theta_hat <- c(beta0_est, beta1_est)
  theta0    <- c(0, 0)
  
  Sigma_hat <- cov(cbind(beta0_boot, beta1_boot))
  if (det(Sigma_hat) < 1e-12) {
    Sigma_inv <- ginv(Sigma_hat)
  } else {
    Sigma_inv <- solve(Sigma_hat)
  }
  
  diff <- theta_hat - theta0
  W    <- as.numeric(t(diff) %*% Sigma_inv %*% diff)
  df   <- 2
  pval <- 1 - pchisq(W, df)
  
  pval
}


#' Wald test for equality of two vaccine-effect parameters.
#'
#' This function tests the null hypothesis
#' \deqn{H_0: \beta_1 = \beta_0}
#' using the bootstrap distribution of the difference
#' \deqn{\delta = \beta_1 - \beta_0}.
#'
#' @param overall_estimates Numeric vector of parameter estimates from the
#'   full dataset (typically from \code{run_corrected_analysis}).
#' @param boot_estimates Numeric matrix whose rows are bootstrap parameter
#'   vectors (typically from \code{run_boot_corrected_analysis}).
#' @param col_beta0 Integer column index in \code{boot_estimates} containing
#'   the \eqn{\beta_0} estimates.
#' @param col_beta1 Integer column index in \code{boot_estimates} containing
#'   the \eqn{\beta_1} estimates.
#'
#' @return A numeric value: the two-sided p-value for the null hypothesis
#'   \eqn{\beta_1 = \beta_0}.
#'
hypothesis_test_equal_VE <- function(overall_estimates, 
                                     boot_estimates, 
                                     col_beta0 = 1, 
                                     col_beta1 = 2){
  
  beta0_est <- overall_estimates[col_beta0]
  beta1_est <- overall_estimates[col_beta1]
  
  beta0_boot <- boot_estimates[,col_beta0]
  beta1_boot <- boot_estimates[,col_beta1]
  
  delta_hat  <- beta1_est - beta0_est
  delta_boot <- beta1_boot - beta0_boot
  se_delta   <- sd(delta_boot)
  
  z_wald     <- delta_hat / se_delta
  p_two_sided <- 2 * pnorm(-abs(z_wald))
  
  p_two_sided
}


#' Wald tests for each vaccine-effect parameter individually.
#'
#' This function performs two separate Wald tests:
#' \deqn{H_{0,0}: \beta_0 = 0}
#' \deqn{H_{0,1}: \beta_1 = 0}
#'
#' The standard errors are estimated from the bootstrap distribution of
#' \eqn{\beta_0} and \eqn{\beta_1}.
#'
#' @param overall_estimates Numeric vector of parameter estimates from the
#'   full data (typically from \code{run_corrected_analysis}).
#' @param boot_estimates Numeric matrix of bootstrap estimates (each row
#'   corresponds to one bootstrap replicate).
#' @param col_beta0 Integer index of the column in \code{boot_estimates}
#'   corresponding to \eqn{\beta_0}.
#' @param col_beta1 Integer index of the column in \code{boot_estimates}
#'   corresponding to \eqn{\beta_1}.
#'
#' @return A list containing:
#' \describe{
#'   \item{beta0_p}{Two-sided p-value for testing \eqn{\beta_0 = 0}.}
#'   \item{beta1_p}{Two-sided p-value for testing \eqn{\beta_1 = 0}.}
#' }
#'
hypothesis_test_VE_0 <- function(overall_estimates, 
                                 boot_estimates, 
                                 col_beta0 = 1, 
                                 col_beta1 = 2) {
  # Extract point estimates
  beta0_est <- overall_estimates[col_beta0]
  beta1_est <- overall_estimates[col_beta1]
  
  # Extract bootstrap replicates
  beta0_boot <- boot_estimates[, col_beta0]
  beta1_boot <- boot_estimates[, col_beta1]
  
  # Compute bootstrap SEs
  se_beta0 <- sd(beta0_boot)
  se_beta1 <- sd(beta1_boot)
  
  # Wald test statistics
  z_beta0 <- beta0_est / se_beta0
  z_beta1 <- beta1_est / se_beta1
  
  # Two-sided p-values
  p_beta0 <- 2 * pnorm(-abs(z_beta0))
  p_beta1 <- 2 * pnorm(-abs(z_beta1))
  
  # Return named vector or list
  out <- list(beta0_p = p_beta0, beta1_p = p_beta1)
  return(out)
}


#' End-to-end analysis: estimation, bootstrap uncertainty, and hypothesis tests.
#'
#' This function runs the full pipeline:
#' \enumerate{
#'   \item Fit the classification-weighted proportional hazards model on the
#'         original data via \code{\link{run_corrected_analysis}}.
#'   \item Obtain bootstrap replicates of the parameter estimates via
#'         \code{\link{run_boot_corrected_analysis}}.
#'   \item Perform three hypothesis tests:
#'         \itemize{
#'           \item \code{\link{hypothesis_test_no_VE}}: joint test of no effect;
#'           \item \code{\link{hypothesis_test_equal_VE}}: equality of the two effects;
#'           \item \code{\link{hypothesis_test_VE_0}}: each effect individually vs 0.
#'         }
#'   \item Compute point estimates and bootstrap-based confidence intervals for
#'         vaccine effectiveness (VE) corresponding to two parameters, indexed
#'         by \code{col_beta0} and \code{col_beta1}.
#' }
#'
#' @param df Data frame containing the analysis variables (event times, status,
#'   covariates, and auxiliary variables).
#' @param cl_probs_vars Character vector of column names used to define groups
#'   for estimating posterior classification probabilities.
#' @param covariates Character vector of covariate names used in the model.
#' @param event_var Name of the event/censoring indicator column.
#' @param time_var Name of the observed time variable.
#' @param stratification_var Optional column name(s) defining strata. If
#'   \code{NULL}, all observations are in a single stratum.
#' @param boots Number of bootstrap samples for uncertainty quantification.
#' @param prior_type Character string passed down to the prior/construction
#'   functions (\code{obtain_prior}, \code{estimate_posterior_probs},
#'   \code{add_classification_probs}).
#' @param degree_spline Spline degree for the prior when
#'   \code{prior_type = "spline"}.
#' @param reg_param Regularization parameter for the spline prior.
#' @param prior_grid Support grid for the nonparametric MLE prior.
#' @param par_init Initial parameter vector passed to \code{run_corrected_analysis}
#'   and ultimately to \code{optim()}.
#' @param thresh Threshold used when computing posterior classification
#'   probabilities.
#' @param col_beta0 Column index (in the parameter vector) corresponding to
#'   the first effect parameter (for VE\eqn{_0}).
#' @param col_beta1 Column index corresponding to the second effect parameter
#'   (for VE\eqn{_1}).
#'
#' @return A list with components:
#' \describe{
#'   \item{ve0_estimate}{Point estimate of VE\eqn{_0 = 1 - \exp(\beta_0)}.}
#'   \item{ve0_se}{Bootstrap standard error for \eqn{\beta_0}.}
#'   \item{ve0_lo}{Lower 95\% confidence limit for VE\eqn{_0}.}
#'   \item{ve0_hi}{Upper 95\% confidence limit for VE\eqn{_0}.}
#'   \item{ve1_estimate}{Point estimate of VE\eqn{_1 = 1 - \exp(\beta_1)}.}
#'   \item{ve1_se}{Bootstrap standard error for \eqn{\beta_1}.}
#'   \item{ve1_lo}{Lower 95\% confidence limit for VE\eqn{_1}.}
#'   \item{ve1_hi}{Upper 95\% confidence limit for VE\eqn{_1}.}
#'   \item{p_ve0}{Two-sided p-value for testing \eqn{\beta_0 = 0}.}
#'   \item{p_ve1}{Two-sided p-value for testing \eqn{\beta_1 = 0}.}
#'   \item{p_val_no_ve}{P-value for the joint no-effect test.}
#'   \item{p_val_equal}{P-value for the equality test \eqn{\beta_1 = \beta_0}.}
#' }
#'
run_full_analysis <- function(df, cl_probs_vars, covariates, 
                              event_var = "allstatus",
                              time_var = "observed_time",
                              stratification_var = NULL,
                              boots = 100,
                              prior_type = "spline",
                              degree_spline = 15, 
                              reg_param = 1, 
                              prior_grid = seq(0, 1, 0.005), 
                              par_init = c(0.69, -0.05), 
                              thresh = 0.01, 
                              col_beta0 = 1, 
                              col_beta1 = 2){
  
  res_overall <- run_corrected_analysis(df, 
                                        cl_probs_vars = cl_probs_vars, 
                                        covariates = covariates, 
                                        event_var = event_var,
                                        time_var = time_var,
                                        stratification_var = stratification_var,
                                        prior_type = prior_type,
                                        degree_spline = degree_spline, 
                                        reg_param = reg_param, 
                                        prior_grid = prior_grid, 
                                        par_init = par_init, 
                                        thresh = thresh)
  
  res_boot <- run_boot_corrected_analysis(df, 
                                          cl_probs_vars = cl_probs_vars, 
                                          covariates = covariates, 
                                          event_var = event_var,
                                          time_var = time_var,
                                          stratification_var = stratification_var,
                                          boots = boots,
                                          prior_type = prior_type,
                                          degree_spline = degree_spline, 
                                          reg_param = reg_param, 
                                          prior_grid = prior_grid, 
                                          par_init = par_init, 
                                          thresh = thresh)
  
  res_boot <- na.omit(res_boot)

  p_val_no_ve <- hypothesis_test_no_VE(res_overall, 
                                          res_boot,
                                         col_beta0 = col_beta0, 
                                         col_beta1 = col_beta1)
  
  p_val_equal <- hypothesis_test_equal_VE(res_overall, 
                                          res_boot,
                                          col_beta0 = col_beta0, 
                                          col_beta1 = col_beta1)
  
  p_vals_0 <- hypothesis_test_VE_0(res_overall, 
                                   res_boot,
                                   col_beta0 = col_beta0, 
                                   col_beta1 = col_beta1)
  
  ve0_estimate <- 1 - exp(res_overall[col_beta0])
  ve0_se <- sd(res_boot[,col_beta0])
  ve0_lo <- 1 - exp(res_overall[col_beta0] - qnorm(0.025) * ve0_se)
  ve0_hi <- 1 - exp(res_overall[col_beta0] + qnorm(0.025) * ve0_se)
  
  ve1_estimate <- 1 - exp(res_overall[col_beta1])
  ve1_se <- sd(res_boot[,col_beta1])
  ve1_lo <- 1 - exp(res_overall[col_beta1] - qnorm(0.025) * ve1_se)
  ve1_hi <- 1 - exp(res_overall[col_beta1] + qnorm(0.025) * ve1_se)

  list(
    ve0_estimate = ve0_estimate, 
    ve0_se = ve0_se, 
    ve0_lo = ve0_lo, 
    ve0_hi = ve0_hi, 
    ve1_estimate = ve1_estimate,
    ve1_se = ve1_se, 
    ve1_lo = ve1_lo,
    ve1_hi = ve1_hi,
    p_ve0 = p_vals_0$beta0_p, 
    p_ve1 = p_vals_0$beta1_p, 
    p_val_no_ve = p_val_no_ve,
    p_val_equal = p_val_equal
  )
}


#' Screen basic variability properties of observed event rates.
#'
#' This function computes simple descriptive measures of the raw event-rate
#' estimates \code{K / sequence_depth}, including the number of observations
#' in each binary category relative to a threshold, the number of observations
#' in each category of \code{mindist_mark}, and how often the raw classification
#' disagrees with \code{mindist_mark}.
#'
#' Specifically, it computes:
#' \itemize{
#'   \item the proportion of observations with raw event rate strictly between 0 and 1;
#'   \item counts of observations with \code{raw_mark = 0} and \code{raw_mark = 1};
#'   \item counts of observations with \code{minddist_mark = 1} and \code{minddist_mark = 0};
#'   \item a mismatch score, the sum of \code{|mindist_mark - (1 - raw_mark)|}.
#' }
#'
#' @param df Data frame containing columns \code{K}, \code{sequence_depth},
#'   and \code{mindist_mark}.
#' @param thresh Numeric threshold for defining the binary variable
#'   \code{raw_mark = I(K / sequence_depth >= thresh)}.
#'
#' @return A list with components:
#' \describe{
#'   \item{prop_other}{Proportion of observations with raw event rate in (0,1).}
#'   \item{cat0_sum}{Count of observations where \code{raw_mark == 0}.}
#'   \item{cat1_sum}{Count of observations where \code{raw_mark == 1}.}
#'   \item{mindist0_sum}{Count of observations where \code{mindist_mark == 1}.}
#'   \item{mindist1_sum}{Count of observations where \code{mindist_mark == 0}.}
#'   \item{diff_from_mindist}{Mismatch score between raw and mindist marks.}
#' }
#'
run_variability_screen <- function(df, thresh = 0.01){
  
  df_new <- df %>% 
    mutate(raw_Q = K / sequence_depth, 
           raw_mark = ifelse(raw_Q >= thresh, 1, 0)) %>% 
    filter(!is.na(raw_Q))
  
  diff_from_mindist <- sum((df_new %>% 
                              mutate(mark_diff = abs(mindist_mark - 1 + raw_mark)))$mark_diff)
  
  cat0_sum = nrow(df_new %>% filter(raw_mark == 0))
  
  cat1_sum = nrow(df_new %>% filter(raw_mark == 1))
  
  list(
    prop_other = sum(df_new$raw_Q > 0 & df_new$raw_Q < 1, 
                     na.rm = T) / nrow(df_new), 
    
    cat0_sum = cat0_sum, 
    cat1_sum = cat1_sum, 
    
    mindist0_sum = sum(df_new$mindist_mark == 1, na.rm =T), 
    mindist1_sum = sum(df_new$mindist_mark == 0, na.rm =T),
    
    diff_from_mindist = diff_from_mindist
  )
  
}


#' Construct an analysis dataset for a specific genetic/sequence mark.
#'
#' This helper function combines the sieve dataset, subject-level dataset,
#' and mindist dataset to produce an analysis-ready data frame containing
#' event times, treatment assignment, covariates, and mark-specific
#' sequencing summary variables (\code{K} and \code{sequence_depth}).
#'
#' @description
#' For the selected \code{mark}, the function:
#' \enumerate{
#'   \item Aggregates sequencing-level data (\code{df_sieve}) to compute, per
#'         subject, the sequencing depth and the count of sequences where the
#'         mark is absent (\code{K = sum(mark == 0)}).
#'   \item Merges these summaries with subject-level data (\code{df_sieve_subj})
#'         containing treatment assignment, event times, and covariates.
#'   \item Sets \code{K} and \code{sequence_depth} to \code{NA} for censored
#'         subjects, defining \code{allstatus = 1} for events and \code{0} for
#'         non-events.
#'   \item Optionally merges mindist-based marks from \code{df_sieve_mindist}
#'         when the selected \code{mark} exists there; otherwise initializes
#'         \code{mindist_mark = 0}.
#' }
#'
#' @param df_sieve A sequencing-level dataset where each row corresponds to a
#'   sequence observation and includes at least \code{subjid} and the
#'   binary mark indicator named by \code{mark}.
#' @param df_sieve_subj A subject-level dataset containing at least:
#'   \code{subjid}, treatment assignment (\code{armdesc}),
#'   event time (\code{hiv1fposday}), event indicator (\code{hiv1event}),
#'   and specified covariates.
#' @param df_sieve_mindist A dataset containing mindist-based binary marks
#'   per subject. Must include \code{subjid} and optionally a column named
#'   exactly \code{mark}.
#' @param mark A character string naming the mark column to use from
#'   \code{df_sieve} (and optionally \code{df_sieve_mindist}).
#'
#' @return A data frame containing:
#' \describe{
#'   \item{subjid}{Subject identifier.}
#'   \item{Z}{Treatment indicator (1 = vaccine, 0 = placebo).}
#'   \item{ind_sa, age, bmi, riskscoresl}{Covariates carried from
#'         \code{df_sieve_subj}.}
#'   \item{observed_time}{Event or censoring time (from \code{hiv1fposday}).}
#'   \item{sequence_depth}{Number of sequences for the selected subject.}
#'   \item{K}{Count of sequences where the mark is absent (mark == 0).}
#'   \item{allstatus}{Event indicator for the analysis:
#'                    1 = event, 0 = censored.}
#'   \item{mindist_mark}{Mindist-based version of the mark, if available,
#'                       otherwise 0.}
#' }
#'
create_mark_dataset <- function(df_sieve, df_sieve_subj, df_sieve_mindist, mark){
  
  df_mark <- df_sieve %>% 
    group_by(subjid) %>% 
    summarise(sequence_depth = n(), 
              K = sum(!!sym(mark) == 0))
  
  df_analysis <- df_sieve_subj %>% 
    
    select(subjid, armdesc, hiv1fposday, hiv1event, ind_sa, age, bmi, riskscoresl) %>% 
    
    merge(df_mark, all.x = T) %>% 
    
    mutate(Z = ifelse(armdesc == "Vaccine", 1, 0), 
           K = ifelse(hiv1event == 1, K, NA),
           sequence_depth = ifelse(hiv1event == 1, sequence_depth, NA),
           allstatus = ifelse(is.na(K), 0, 1)) %>% 
    
    # excludes missing marks
    filter(!(hiv1event == 1 & allstatus == 0)) %>% 
    
    rename(observed_time = hiv1fposday) %>% 
    
    select(subjid, Z, ind_sa, age, bmi, riskscoresl, 
           observed_time, sequence_depth, K, allstatus)
  
  if (mark %in% names(df_sieve_mindist)) {
    df_analysis <- df_analysis %>% 
      merge(df_sieve_mindist %>% 
              select(subjid, !!sym(mark)), 
            by = "subjid", 
            all.x = T) %>% 
      rename(mindist_mark = !!sym(mark))
  } else {
    df_analysis <- df_analysis %>% 
      mutate(mindist_mark = 0)
  }
  
  df_analysis
}