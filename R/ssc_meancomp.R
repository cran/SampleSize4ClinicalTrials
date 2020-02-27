#' @title Sample Size Calculation for the Comparison of Means in Phase III Clinical Trials
#'
#' @description This function aims to calculate sample size for the comparison of means
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
#' @param sigma The variance of observed outcomes in both arms
#' @param theta The difference between means in the two arms
#' @param delta The prespecified non-inferiority or equivalence margin in non-inferiority (3) or equivalence (4) trials
#'
#' @return samplesize
#'
#' @usage ssc_meancomp(design = c(1,2,3,4), ratio, alpha, power, sigma, theta, delta)
#'
#' @importFrom stats qnorm
#'
#' @examples
#' ##The comparison of means, a non-inferiority trial with a non-inferiority margin 0.5
#' ##the true treatment difference is assumed to be zero in non-inferiority and equivalence trials
#' ssc_meancomp(design=3, ratio=1, alpha=0.05, power=0.9, sigma=1, theta=0, delta=0.5)
#'
#' @references Yin, G. (2012). Clinical Trial Design: Bayesian and Frequentist Adaptive Methods. John Wiley & Sons.
#'
#' @export

##Sample size calculation for mean comparison
ssc_meancomp<-function(design = c(1,2,3,4), ratio, alpha, power, sigma ,theta, delta){
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

  ##Sample size equation, n2 means number of subjects in the control arm
  n2<-ceiling((1+1/ratio)*(sigma^2)*(z_alpha+z_beta)^2/denom)
  n1<-ratio*n2

  ##Calculate the sample size
  samplesize<- data.frame(n1, n2)
  colnames(samplesize) <- c("Treatment", "Control")
  return(samplesize)
}
