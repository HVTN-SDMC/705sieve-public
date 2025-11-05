##################################################
### Calculate the adjusted p-values to control familywise error rate(FWER)and false discovery rate(FDR)
### using the resampleing based methods.
###
### Input:
###     p.unadj: an 1xm vector of unjected p-values calculated from the original data set
###     p.perms:  an Bxm matrix of p-values calculated from B sets of data sets that are resampled from the orginal data
###         set under null hypothese
###   alpha:  any unjected p-values less than alpha will not be calculated for adjusted p-values and their adjusted
###           p-values are NA.
### Output:
### p.FWER: an 1xm vector of adjusted p-values to control FWER
### p.FDR:  an 1xm vector of adjusted p-values to control FDR
###
### FWER adjusted p-values, P.FWER, are calculated based on the resampling step down procedure (Westfall and Young 1993).
### FDR adjusted p-values, p.FDR, are calculated based on the estimations of E(R0)/E(R) where E(R0) is the
###     expectation of the number of rejected null hypotheses (R0) that is estimated from the resampled data sets
###     under the null hypotheses; E(R) is the expection of the number of all rejected hypotheses (R) that is estimated
###     by the maximum of R0 and the number of rejections from the observed data set.
### According to Jensen inequality and R0 is a positive linear function of R,
###   FDR=E(R0/R) <= E(R0)/E(R). Therefore, our estimation for p.FDR would control below the FDR level.
### Ref:
###     Westfall and Young 1993 "Resampling-based multiple testing: Examples and methods for p-value
###         adjustment", John Wiley & Sons.
###     Westfall and Troendle "Multiple testing wirh minimum aasumptions", Biom J. 2008
###     Storey and Tibshrani "Statistical siggnificance for genomewide studies", PNAS 2003
###
### Created by Sue Li, 4/2015
### Modified by Sue Li, 5/31/2023
##################################################
p.adj.perm2 <- function(p.unadj,p.perms,alpha=0.05)
{
  stopifnot( ncol( p.perms ) == length( p.unadj ) );
  
  # If "p.unadj" has no names, give it names either from colnames( p.perms ) or as 1:length( p.unadj ).
  if( is.null( names( p.unadj ) ) ) {
    if( is.null( colnames( p.perms ) ) ) {
      names( p.unadj ) <- 1:length( p.unadj );
    } else {
      names( p.unadj ) <- colnames( p.perms );
    }
  }
  B = dim(p.perms)[1]
  m = length(p.unadj)

  ### order p from the smallest to the largest
  # We must sort the p.unadj values first.
  mode( p.unadj ) <- "numeric";
  which.are.NA <- which( is.na( p.unadj ) );
  p.unadj.order <- order( p.unadj ); # Note this puts NAs last/"largest".
  # We must maintain the same ordering between p.unadj and the columns of p.perms.
  p.unadj <- p.unadj[ p.unadj.order ];
  p.perms <- p.perms[ , p.unadj.order, drop = FALSE ];
  
  len = sum(round(p.unadj,2)<=alpha)
  # calculate FWER-adjusted p-values
  p.FWER=rep(NA,length(p.unadj))
  
  for (j in 1:len)
  {
      p.FWER[j] = sum((apply(p.perms[,j:m,drop=FALSE], 1, min, na.rm = T)<=p.unadj[j]))/B
  }
  ## enforce monotonicity using successive maximization
  p.FWER[1:len] = cummax(p.FWER[1:len])
  
  # calculate FDR-adjusted p-values
  p.FDR=rep(NA,length(p.unadj))
  for (j in 1:len)
  {
      ## given each p-value
      ## estimate the expectation of # of rejections under null hypotheses
      R0_by_resample = apply(p.perms<=p.unadj[j], 1, sum, na.rm = T )
      ER0 = sum(R0_by_resample)/B
      ## calculate # of rejections observed in the data
      # R.ob = j  #sum(p.unadj<=p.unadj[j])
      R.ob = sum(p.unadj<=p.unadj[j])
      ## R is max(R0,R)
      R = sum(pmax(R0_by_resample,R.ob))/B
      ## FDR=E(R0/R|R>0) FDR=0 if R=0
      p.FDR[j]=min(ifelse(R>0,ER0/R,0),1)
  }
  
  o1=order(p.unadj[1:len],decreasing=TRUE)
  ro=order(o1)
  
  p.FDR[1:len]=pmin(1,cummin(p.FDR[1:len][o1]))[ro]
  
  ## the results are in an ascending order of the unadjusted p-values
  .rv <- cbind(p.unadj,p.FWER,p.FDR)
  rownames(.rv) <- names(p.unadj)
  return(.rv)
}
  
