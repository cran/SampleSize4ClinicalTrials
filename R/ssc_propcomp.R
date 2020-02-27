#' @title Sample Size Calculation for the Comparison of Proportions in Phase III Clinical Trials
#'
#' @description This function aims to calculate sample size for the comparison of proportions
#' in Phase III clinical trials.
#'
#' @param design The design of the clinical trials.
#'   \cr 1
#'     \cr Testing for equality
#'   \cr 2
#'     \cr Superiority trial
#'   \cr 3
#'     \cr Non-inferiority trial
#'   \cr 4
#'     \cr Equivalence trial.
#' @param ratio The ratio between the number of subjects in the treatment arm and that in the control arm.
#' @param alpha Type I error rate
#' @param power Statistical power of the test (1-type II error rate)
#' @param p1 The response rate of the treatment arm
#' @param p2 The response rate of the control arm
#' @param theta The difference between proportions in the two arms
#' @param delta The prespecified non-inferiority or equivalence margin in non-inferiority (3) or equivalence (4) trials
#'
#' @return samplesize
#'
#' @usage ssc_propcomp(design=c(1,2,3,4), ratio, alpha, power, p1, p2, theta, delta)
#'
#' @importFrom stats qnorm
#'
#' @examples
#' ##The comparison of proportions, a superiority trial
#' ssc_propcomp(design=2, ratio=3, alpha=0.025, power=0.8, p1=0.4, p2=0.2, theta=0.2)
#' @references Yin, G. (2012). Clinical Trial Design: Bayesian and Frequentist Adaptive Methods. John Wiley & Sons.
#'
#' @export

##Sample size calculation for the comparison of proportions
ssc_propcomp<-function(design = c(1,2,3,4), ratio, alpha, power, p1, p2, theta, delta){
  if (!(design %in% 1:4))
    stop("Unrecognized study design, please select from: 1. Test for equility, 2. Superiority trial,
         3. Non-inferiority trial, 4. Equivalence trial")
  ##Assign values to z_alpha
  if (design==1)
    z_alpha<-qnorm(1-(alpha/2))
  else
    z_alpha<-qnorm(1-alpha)
  ##Assign values to z_beta
  if (design == 4)
    z_beta<-qnorm((1+power)/2)
  else
    z_beta<-qnorm(power)
  ##Denominator
  if (design %in% 1:2)
    denom<-theta^2
  if (design==3)
    denom<-(theta+delta)^2
  if (design==4)
    denom<-(delta-abs(theta))^2

  ##Pooled variance, n4 means number of sujects in the control arm
  pbar<-(p1+p2)/2
  n4<-ceiling(((1+1/ratio)*(z_alpha*sqrt(pbar*(1-pbar))+z_beta*sqrt((p1*(1-p1)+p2*(1-p2)*ratio)/(ratio+1)))^2)/denom)
  n3<-ratio*n4

  ##Calculate the sample size
  samplesize<- data.frame(n3, n4)
  colnames(samplesize) <- c("Treatment", "Control")
  return(samplesize)
}
