input_check <- function(Y,X,B,A,a,b){
  if(!is.numeric(A)||!is.vector(A)){
    stop("Please enter a numeric vector for treatment assignments (A)!")
  }
  if(!is.numeric(B)||!is.vector(B)){
    stop("Please enter a numeric vector for stratum labels (B)!")
  }
  if(!is.numeric(Y)||!is.vector(Y)){
    stop("Please enter a numeric vector for observed outcomes (Y)!")
  }
 if(is.null(X)){
    if(length(A)!=length(B)||length(A)!=length(Y)){
      stop("A, B, Y should have same lengths!")
    }
  }
  else if(!is.numeric(X)||!is.matrix(X)){
    stop("Please enter a numeric matrix for additional covariates!")
  
    if(length(A)!=length(B)||length(A)!=length(Y)||length(A)!=nrow(X)){
      stop("Lengths of A, B, Y and the number of rows in X should be the same!")
    }
  }
  
  if((!(a %in% A)) || (!(b %in% A))){
    stop("Treatment a and b should be in the treatment assignment A!")
  }
  if(a==b){
    stop("Please enter two different treatments!")
  }
}

#' Plug-in estimator
#' 
#' Estimating and inferring the treatment effect based on the plug-in estimator.
#' 
#' Estimating and inferring the treatment effect based on the plug-in estimator. The estimator can also be
#' derived by fully saturated model in Bugni et al. (2019). It is the benchmark estimator described in Section 3.1, 
#' Gu et al. (2023).
#' 
#' @param Y a numeric vector of observed outcomes. Its length should be the same as the number of subjects.
#' @param A a numeric vector of treatment assignments. The assignments should be consecutive integers from 0 to nA, where 
#' nA is the number of different treatments (including control group). 
#' @param B a numeric vector of stratum labels. Its length should be the same as the number of subjects.
#' @param a,b an integer denote the treatment in \code{A}. The result will present the estimated treatment effect between selected
#' treatment \code{a} and treatment \code{b}
#' @param conf.level confidence level of the interval. Default is 0.95.
#' 
#' @return A list of class \code{"htest"} containing the following components:
#' \item{statistic}{the value of the t-statistic.} \item{p.value}{the p-value for
#' the test.} \item{conf.int}{a confidence interval under chosen level
#' \code{conf.level} for the difference in treatment effect between treatment
#' group and control group.} \item{estimate}{estimated treatment effect
#' difference between treatment group and control group.} \item{method}{a
#' character string indicating what type of regression was performed.}
#' 
#' @references Bugni, F. A., Canay, I. A., & Shaikh, A. M. (2019). \emph{Inference under covariate-adaptive 
#' randomization with multiple treatments}. Quantitative Economics 10, 1747–1785.
#' @references Gu Y., Liu H., Ma W.(2023). \emph{Regression-based multiple treatment effect estimation under covariate-adaptive 
#' randomization.} In review.
#' 
#' @examples 
#' #The code replicates the simulation setting of Model 1 in Section 4, Gu et al. (2023).
#' 
#' nA = 3
#' nS = 4
#' N=600
#' p_s = c(2,2,3,3)/10
#' 
#' pi_a = rep(c(1,1,1)/3)
#' 
#' p=1
#' beta0 <- matrix(rep(1,p),p,1)
#' beta1 <- matrix(rep(2,p),p,1)
#' beta2 <- matrix(rep(3,p),p,1)
#' 
#' X1 <- sample(1:nS,N,prob = p_s,replace = TRUE)
#' #additional variable
#' X2 <- matrix(rnorm(N*p,0,1),N,p)
#' 
#' trt <- sample(0:(nA-1),N,replace=TRUE,pi_a)
#' 
#' ind0 <- which(trt==0)
#' ind1 <- which(trt==1)
#' ind2 <- which(trt==2)
#' 
#' 
#' Y0 <- 0 + X1[ind0] + X2[ind0,]%*%beta0
#' Y1 <- 1 + X1[ind1] + X2[ind1,]%*%beta1
#' Y2 <- 2 + X1[ind2] + X2[ind2,]%*%beta2
#' Y <- rep(0,N)
#' Y[ind0] <- Y0
#' Y[ind1] <- Y1
#' Y[ind2] <- Y2
#' Y <- Y + rnorm(N,0,1)
#' tau.pl(Y,A=trt,B=X1,a=1,b=0)
#' @export
tau.pl <- function(Y,A,B,a,b,conf.level=0.95){
  input_check(Y,X=NULL,B,A,a,b)
  N <- length(Y)
  
  tau_vec <- tau_sat(Y,B,A)
  estimate <- tau_vec[a+1] - tau_vec[b+1]
  var_pl <- sd_sat(Y,B,A,a,b)
  
  stderr <- sqrt(var_pl/N)
  testmethod <- "Benchmark estimator"
  tstat <- estimate/stderr
  pval <- dplyr::if_else(stats::pnorm(tstat)<(1-stats::pnorm(tstat)),stats::pnorm(tstat),1-stats::pnorm(tstat))*2
  cint <- c(estimate + stderr*stats::qnorm((1-conf.level)/2),estimate - stderr*stats::qnorm((1-conf.level)/2))
  attr(cint,"conf.level") <- conf.level
  names(tstat) <- "t"
  names(estimate) <- "difference in treatment effect"
  rval<-list(statistic = tstat, p.value = pval, conf.int = cint,
             estimate = estimate, stderr = stderr,
             method = testmethod)
  class(rval) <- "htest"
  return(rval)
}


#' Stratum-common estimator
#' 
#' Estimating and inferring the treatment effect based on the stratum-common estimator.
#' 
#' Estimating and inferring the treatment effect based on the stratum-common estimator. It implements the 
#' methods as described in Section 3.2, Gu et al.(2023).
#' 
#' @param Y a numeric vector of observed outcomes. Its length should be the same as the number of subjects.
#' @param X a numeric design matrix containing additional covariates used in the model.
#' @param A a numeric vector of treatment assignments. The assignments should be consecutive integers from 0 to nA, where 
#' nA is the number of different treatments (including control group). 
#' @param B a numeric vector of stratum labels. Its length should be the same as the number of subjects.
#' @param a,b an integer denote the treatment in \code{A}. The result will present the estimated treatment effect between selected
#' treatment \code{a} and treatment \code{b}
#' @param conf.level confidence level of the interval. Default is 0.95.
#' 
#' @return A list of class \code{"htest"} containing the following components:
#' \item{statistic}{the value of the t-statistic.} \item{p.value}{the p-value for
#' the test.} \item{conf.int}{a confidence interval under chosen level
#' \code{conf.level} for the difference in treatment effect between treatment
#' group and control group.} \item{estimate}{estimated treatment effect
#' difference between treatment group and control group.} \item{method}{a
#' character string indicating what type of regression was performed.}
#' 
#' @references Gu Y., Liu H., Ma W.(2023). \emph{Regression-based multiple treatment effect estimation under covariate-adaptive 
#' randomization.} In review.
#' 
#' @examples 
#' #The code replicates the simulation setting of Model 1 in Section 4, Gu et al. (2023).
#' 
#' nA = 3
#' nS = 4
#' N=600
#' p_s = c(2,2,3,3)/10
#' 
#' pi_a = rep(c(1,1,1)/3)
#' 
#' p=1
#' beta0 <- matrix(rep(1,p),p,1)
#' beta1 <- matrix(rep(2,p),p,1)
#' beta2 <- matrix(rep(3,p),p,1)
#' 
#' X1 <- sample(1:nS,N,prob = p_s,replace = TRUE)
#' #additional variable
#' X2 <- matrix(rnorm(N*p,0,1),N,p)
#' 
#' trt <- sample(0:(nA-1),N,replace=TRUE,pi_a)
#' 
#' ind0 <- which(trt==0)
#' ind1 <- which(trt==1)
#' ind2 <- which(trt==2)
#' 
#' 
#' Y0 <- 0 + X1[ind0] + X2[ind0,]%*%beta0
#' Y1 <- 1 + X1[ind1] + X2[ind1,]%*%beta1
#' Y2 <- 2 + X1[ind2] + X2[ind2,]%*%beta2
#' Y <- rep(0,N)
#' Y[ind0] <- Y0
#' Y[ind1] <- Y1
#' Y[ind2] <- Y2
#' Y <- Y + rnorm(N,0,1)
#' tau.sc(Y,X2,trt,X1,1,0)
#' @export
tau.sc <- function(Y,X,A,B,a,b,conf.level=0.95){
  input_check(Y,X,B,A,a,b)
  N <- length(Y)
  
  tau_vec <- tau_sc(Y,X,B,A)
  estimate <- tau_vec[a+1] - tau_vec[b+1]
  var_sc <- sd_sc_ne(Y,X,B,A,a,b)
  stderr <- sqrt(var_sc/N)
  testmethod <- "Stratum common estimator"
  tstat <- estimate/stderr
  pval <- dplyr::if_else(stats::pnorm(tstat)<(1-stats::pnorm(tstat)),stats::pnorm(tstat),1-stats::pnorm(tstat))*2
  cint <- c(estimate + stderr*stats::qnorm((1-conf.level)/2),estimate - stderr*stats::qnorm((1-conf.level)/2))
  attr(cint,"conf.level") <- conf.level
  names(tstat) <- "t"
  names(estimate) <- "difference in treatment effect"
  rval<-list(statistic = tstat, p.value = pval, conf.int = cint,
             estimate = estimate, stderr = stderr,
             method = testmethod)
  class(rval) <- "htest"
  return(rval)
}


#' Stratum-specific estimator
#' 
#' Estimating and inferring the treatment effect based the on the stratum-specific estimator.
#' 
#' Estimating and inferring the treatment effect based on the stratum-specific estimator. It implements the 
#' methods as described in Section 3.3, Gu et al.(2023).
#' 
#' @param Y a numeric vector of observed outcomes. Its length should be the same as the number of subjects.
#' @param X a numeric design matrix containing additional covariates used in the model.
#' @param A a numeric vector of treatment assignments. The assignments should be consecutive integers from 0 to nA, where 
#' nA is the number of different treatments (including control group). 
#' @param B a numeric vector of stratum labels. Its length should be the same as the number of subjects.
#' @param a,b an integer denote the treatment in \code{A}. The result will present the estimated treatment effect between selected
#' treatment \code{a} and treatment \code{b}
#' @param conf.level confidence level of the interval. Default is 0.95.
#' 
#' @return A list of class \code{"htest"} containing the following components:
#' \item{statistic}{the value of the t-statistic.} \item{p.value}{the p-value for
#' the test.} \item{conf.int}{a confidence interval under chosen level
#' \code{conf.level} for the difference in treatment effect between treatment
#' group and control group.} \item{estimate}{estimated treatment effect
#' difference between treatment group and control group.} \item{method}{a
#' character string indicating what type of regression was performed.}
#' 
#' @references Gu Y., Liu H., Ma W.(2023). \emph{Regression-based multiple treatment effect estimation under covariate-adaptive 
#' randomization.} In review.
#' 
#' @examples 
#' #The code replicates the simulation setting of Model 1 in Section 4, Gu et al. (2023).
#' 
#' nA = 3
#' nS = 4
#' N=600
#' p_s = c(2,2,3,3)/10
#' 
#' pi_a = rep(c(1,1,1)/3)
#' 
#' p=1
#' beta0 <- matrix(rep(1,p),p,1)
#' beta1 <- matrix(rep(2,p),p,1)
#' beta2 <- matrix(rep(3,p),p,1)
#' 
#' X1 <- sample(1:nS,N,prob = p_s,replace = TRUE)
#' #additional variable
#' X2 <- matrix(rnorm(N*p,0,1),N,p)
#' 
#' trt <- sample(0:(nA-1),N,replace=TRUE,pi_a)
#' 
#' ind0 <- which(trt==0)
#' ind1 <- which(trt==1)
#' ind2 <- which(trt==2)
#' 
#' 
#' Y0 <- 0 + X1[ind0] + X2[ind0,]%*%beta0
#' Y1 <- 1 + X1[ind1] + X2[ind1,]%*%beta1
#' Y2 <- 2 + X1[ind2] + X2[ind2,]%*%beta2
#' Y <- rep(0,N)
#' Y[ind0] <- Y0
#' Y[ind1] <- Y1
#' Y[ind2] <- Y2
#' Y <- Y + rnorm(N,0,1)
#' tau.ss(Y,X2,trt,X1,1,0)
#' @export
tau.ss <- function(Y,X,A,B,a,b,conf.level=0.95){
  input_check(Y,X,B,A,a,b)
  N <- length(Y)
  
  tau_vec <- tau_ss(Y,X,B,A)
  estimate <- tau_vec[a+1] - tau_vec[b+1]
  var_ss <- sd_ss(Y,X,B,A,a,b)
  
  stderr <- sqrt(var_ss/N)
  testmethod <- "Stratum specific estimator"
  tstat <- estimate/stderr
  pval <- dplyr::if_else(stats::pnorm(tstat)<(1-stats::pnorm(tstat)),stats::pnorm(tstat),1-stats::pnorm(tstat))*2
  cint <- c(estimate + stderr*stats::qnorm((1-conf.level)/2),estimate - stderr*stats::qnorm((1-conf.level)/2))
  attr(cint,"conf.level") <- conf.level
  names(tstat) <- "t"
  names(estimate) <- "difference in treatment effect"
  rval<-list(statistic = tstat, p.value = pval, conf.int = cint,
             estimate = estimate, stderr = stderr,
             method = testmethod)
  class(rval) <- "htest"
  return(rval)
}