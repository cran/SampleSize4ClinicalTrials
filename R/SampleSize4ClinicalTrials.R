#' @title Sample Size Calculation for Mean and Proportion Comparisons in Phase 3 Clinical Trials
#'
#' @description This package aims to calculate sample size for mean and proportion
#' comparisons in phase 3 clinical trials.
#'
#' @param cat Type of the outcome for comparison.
#'  \cr "m"
#'    \cr stands for means comparison.
#'  \cr "p"
#'    \cr stands for proportions comparison.
#' @param design The design of the clinical trials.
#'   \cr 1
#'     \cr Testing for equality
#'   \cr 2
#'     \cr Superiority trial
#'   \cr 3
#'     \cr Non-inferiority trial
#'   \cr 4
#'     \cr Equivalence trial.
#' @param ratio The ratio between number of subjects in experimental drug and that in standard drug.
#' @param alpha Type 1 error rate
#' @param power Statistical power of the test (1-type 2 error rate)
#' @param sigma The variance of observed outcomes in both arms (specified in means comparison for continuous outcomes)
#' @param p1 The response rate of the experimental drug (specified proportion comparison for binary outcomes)
#' @param p2 The response rate of the standard drug
#' @param theta The difference between means or proportions in two arms
#' @param delta The prespecified non-inferiority or equivalence margin in non-inferiority or equivalence trials
#'
#' @return samplesize
#'
#' @usage ssc(cat=c("m","p"),design=c(1,2,3,4),ratio,alpha,power,sigma,p1,p2,theta,delta)
#'
#' @importFrom stats qnorm
#'
#' @examples
#' ##Means comparison for continuous outcomes, a non-inferiority trial with a non-inferiority margin 0.5
#' ##the true treatment difference is assumed to be zero in non-inferiority and equivalence trials
#' ssc(cat="m",design=3,ratio=1,alpha=0.05,power=0.9,sigma=1,theta=0,delta=0.5)
#' ##Proportions comparison for binary outcomes, a superiority trial
#' ssc(cat="p",design=2,ratio=3,alpha=0.025,power=0.8,p1=0.4,p2=0.2,theta=0.2)
#' @references Yin, G. (2012). Clinical Trial Design: Bayesian and Frequentist Adaptive Methods. John Wiley & Sons.
#'
#' @export

###sample size calculation
ssc<-function(cat=c("m","p"),design=c(1,2,3,4),ratio,alpha,power,sigma=0,p1=0,p2=0,theta,delta){
  ##assign values to zalpha
  if (design==1)
    zalpha<-qnorm(1-(alpha/2))
  else
    zalpha<-qnorm(1-alpha)
  ##assign values to zbeta
  if (design<4)
    zbeta<-qnorm(power)
  else
    zbeta<-qnorm((1+power)/2)
  ##denominator
  if (design<3)
    denom<-theta^2
  if (design==3)
    denom<-(theta+delta)^2
  if (design==4)
    denom<-(delta-abs(theta))^2
  ##sample size equation, n2 means number of subjects in standard treatment group
  n2<-ceiling((1+1/ratio)*(sigma^2)*(zalpha+zbeta)^2/denom)
  n1<-ratio*n2
  pbar<-(p1+p2)/2
  ##pooled variance, n4 means number of sujects in standard treatment group
  n4<-ceiling(((1+1/ratio)*(zalpha*sqrt(pbar*(1-pbar))+zbeta*sqrt((p1*(1-p1)+p2*(1-p2)*ratio)/(ratio+1)))^2)/denom)
  n3<-ratio*n4
  if (cat=="m")
    samplesize<-list(Category="Sample size for means comparison",Experimental=c("Experimental drug:",n1),Standard=c("Standard drug:",n2))
  else
    samplesize<-list(Category="Sample size for proportions comparison",Experimental=c("Experimental drug:",n3),Standard=c("Standard drug:",n4))
  return(samplesize)
}






