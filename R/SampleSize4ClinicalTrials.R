#' @title Sample Size Calculation for the Comparison of Means or Proportions in Phase III Clinical Trials
#'
#' @description This function aims to calculate sample size for the comparison of means or proportions
#' in Phase III clinical trials.
#'
#' @param cat Type of the outcome for comparison.
#'  \cr "m"
#'    \cr stands for the comparison of means.
#'  \cr "p"
#'    \cr stands for the comparison of proportions.
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
#' @param sigma The variance of observed outcomes in both arms (specified in the comparison of means for continuous outcomes)
#' @param p1 The response rate of the treatment arm (specified the comparison of proportions for binary outcomes)
#' @param p2 The response rate of the control arm
#' @param theta The difference between means or proportions in the two arms
#' @param delta The prespecified non-inferiority or equivalence margin in non-inferiority or equivalence trials
#'
#' @return samplesize
#'
#' @usage ssc(cat=c("m","p"), design=c(1,2,3,4), ratio, alpha, power, sigma, p1, p2, theta, delta)
#'
#' @importFrom stats qnorm
#'
#' @examples
#' ##The comparison of means, a non-inferiority trial with a non-inferiority margin 0.5
#' ##the true treatment difference is assumed to be zero in non-inferiority and equivalence trials
#' ssc(cat="m",design=3,ratio=1,alpha=0.05,power=0.9,sigma=1,theta=0,delta=0.5)
#' ##The comparison of proportions, a superiority trial
#' ssc(cat="p",design=2,ratio=3,alpha=0.025,power=0.8,p1=0.4,p2=0.2,theta=0.2)
#' @references Yin, G. (2012). Clinical Trial Design: Bayesian and Frequentist Adaptive Methods. John Wiley & Sons.
#'
#' @export

##Sample size calculation
ssc<-function(cat=c("m","p"),design=c(1,2,3,4),ratio,alpha,power,sigma=0,p1=0,p2=0,theta,delta){
  ##Assign values to z_alpha
  if (design==1)
    z_alpha<-qnorm(1-(alpha/2))
  else
    z_alpha<-qnorm(1-alpha)
  ##Assign values to z_beta
  if (design<4)
    z_beta<-qnorm(power)
  else
    z_beta<-qnorm((1+power)/2)
  ##Denominator
  if (design<3)
    denom<-theta^2
  if (design==3)
    denom<-(theta+delta)^2
  if (design==4)
    denom<-(delta-abs(theta))^2
  ##Sample size equation, n2 means number of subjects in the control arm
  n2<-ceiling((1+1/ratio)*(sigma^2)*(z_alpha+z_beta)^2/denom)
  n1<-ratio*n2
  pbar<-(p1+p2)/2
  ##Pooled variance, n4 means number of sujects in the control arm
  n4<-ceiling(((1+1/ratio)*(z_alpha*sqrt(pbar*(1-pbar))+z_beta*sqrt((p1*(1-p1)+p2*(1-p2)*ratio)/(ratio+1)))^2)/denom)
  n3<-ratio*n4
  if (cat=="m")
    samplesize<-list(Category="Sample size for the comparison of means",Treatment=c("The treatment arm:",n1),Control=c("The control arm:",n2))
  else
    samplesize<-list(Category="Sample size for the comparison of comparisons",Treatment=c("The treatment arm:",n3),Control=c("The control arm:",n4))
  return(samplesize)
}
