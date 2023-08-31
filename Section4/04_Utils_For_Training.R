#-----------------------------------------------------------------------------|
# SECTION 4.2.4 Simulating NBA Plus-Minus scores without changepoint          |
#                                                                             |
#-----------------------------------------------------------------------------|

################################################################################
##  Functions for Simulating the Plus-Minus score of matches                  ##
##  Pre-processing 1:                                                         ##
##  Modeling Dependencies and functions for Time Series generation            ##
################################################################################

# Description--------------------------------------------------------|
#
#-------------------------------------------------------------------------|

# Remark------------------------------------------------------------------|
#                                                                         |
#-------------------------------------------------------------------------|

#------------------------------------------------------------------------|
# Modeling Dependencies--------------------------------------------------|
add_fn <- list()

# Without dependence
add_fn[["None"]]  <- numeric

# First scenario: AR(1) with coefficient 0.6 (x 1, 2, 4 and 8).
add_fn[["AR1x1"]]  <- function(n) 1*arima.sim(n=n, list(ar=0.6, sd=0.5))
add_fn[["AR1x2"]]  <- function(n) 2*arima.sim(n=n, list(ar=0.6, sd=0.5))
add_fn[["AR1x4"]]  <- function(n) 4*arima.sim(n=n, list(ar=0.6, sd=0.5))
add_fn[["AR1x8"]]  <- function(n) 8*arima.sim(n=n, list(ar=0.6, sd=0.5))

# Second scenario: sin waved with a period of 20 matches (x 1, 2, 4 or 8).
add_fn[["SINx1"]]   <- function(n)  1*sin(seq(0, n/20, length.out=n)*2*pi)
add_fn[["SINx2"]]   <- function(n)  2*sin(seq(0, n/20, length.out=n)*2*pi)
add_fn[["SINx4"]]   <- function(n)  4*sin(seq(0, n/20, length.out=n)*2*pi)
add_fn[["SINx8"]]   <- function(n)  8*sin(seq(0, n/20, length.out=n)*2*pi)


#functions---------------------------------------------------------------|
# Time Series Modeling (with Dependencies)-------------------------------|
sim_fn <- list()
#------------------------------------------------------------------------|
#' @title sim_fn[["Skellam"]]
#'
#' @description Skellam distribution :  mean 1 + dependence
#' @param n number of data points
#' @param mu_skellam mean parameter
#' @param add_dep type of the dependence (AR or sin)
#' @return time series
#'
sim_fn[["Skellam"]] <- function(n, mu_skellam, add_dep) {
  add_to_mu1 <- add_dep(n)
  rskellam(n, mu1 = mu_skellam + add_to_mu1, mu2 = mu_skellam)
}

#------------------------------------------------------------------------|
#' @title sim_fn[["NegBin-500"]]
#'
#' @description Difference of two Negative Binomial distributions with scale parameter neg.size = 500 : mean 1 + dependence
#' @param n number of data points
#' @param mu_skellam mean parameter
#' @param add_dep type of the dependence (AR or sin)
#' @return time series
#'
sim_fn[["NegBin-500"]] <- function(n, mu_skellam, add_dep, neg.size = 500) {
  add_to_mu1 <- add_dep(n)
  rnbinom(n = n, size = neg.size, mu = mu_skellam + add_to_mu1) - rnbinom(n = n, size = neg.size, mu = mu_skellam)
}

#------------------------------------------------------------------------|
#' @title sim_fn[["NegBin-100"]]
#'
#' @description Difference of two Negative Binomial distributions with scale parameter neg.size = 100 :  mean 1 + dependence
#' @param n number of data points
#' @param mu_skellam mean parameters
#' @param add_dep type of the dependence (AR or sin)
#' @return time series
#'
sim_fn[["NegBin-100"]] <- function(n, mu_skellam, add_dep, neg.size = 100) {
  add_to_mu1 <- add_dep(n)
  rnbinom(n = n, size = neg.size, mu = mu_skellam + add_to_mu1) - rnbinom(n = n, size = neg.size, mu = mu_skellam)
}

#------------------------------------------------------------------------|
#' @title sim_fn[["NegBin-50"]]
#'
#' @description Difference of two Negative Binomial distributions with scale parameter neg.size = 50 : mean 1 + dependence
#' @param n number of data points
#' @param mu_skellam mean parameters
#' @param add_dep type of the dependence (AR or sin)
#' @return time series
#'
sim_fn[["NegBin-50"]] <- function(n, mu_skellam, add_dep, neg.size = 50) {
  add_to_mu1 <- add_dep(n)
  rnbinom(n = n, size = neg.size, mu = mu_skellam + add_to_mu1) - rnbinom(n = n, size = neg.size, mu = mu_skellam)
}

#------------------------------------------------------------------------|
#' @title sim_fn[["NegBin-10"]]
#'
#' @description Difference of Negative Binomial distributions with scale parameter neg.size = 10 : mean 1 + dependence
#' @param n number of data points
#' @param mu_skellam mean parameters
#' @param add_dep type of the dependence (AR or sin)
#' @return time series
#'
sim_fn[["NegBin-10"]] <- function(n, mu_skellam, add_dep, neg.size = 10) {
  add_to_mu1 <- add_dep(n)
  rnbinom(n = n, size = neg.size, mu = mu_skellam + add_to_mu1) - rnbinom(n = n, size = neg.size, mu = mu_skellam)
}

#useful functions-------------------------------------------------------|
#-----------------------------------------------------------------------|
#' @title getMeanSd
#'
#' @description
#' @param reg
#' @param x
#' @return
getMeanSd <- function(reg, x){
  n <- length(x)
  out <- c(mean(x[reg]), mean(x[reg] > 0), sd(x[reg]))
  names(out) <- c("Mean PTS Diff", "Mean Win Rate", "Std. Dev")
  out
}
#------------------------------------------------------------------------|
#' @title computeLikeSeg
#'
#' @description
#' @param reg
#' @param x
#' @return
#'
computeLikeSeg <- function(x, reg){
  length(reg)*log(sum( ( x[reg] - mean(x[reg]) )^2 ) / length(reg))
}
#------------------------------------------------------------------------|
#' @title computeLikeAll
#'
#' @description
#' @param reg
#' @param x
#' @return
#'
computeLikeAll <- function(tau, x){
  n <- length(x)
  computeLikeSeg(x, 1:tau) + computeLikeSeg(x, (tau+1):n) - computeLikeSeg(x, 1:n)
}

#------------------------------------------------------------------------|
#' @title sim_Max_ChangeInMean
#'
#' @description Function to simulate 1d data and to return the max statistics of MdFOCuS for Gaussian Model(change in mean)
#' @param n number of data points
#' @param sim.data distribution (by default, rnorm)
#' @return max statistics

sim_Max_ChangeInMean <- function(n, sim.data = rnorm) {
 data <- matrix(sim.data(n), nc = 1) ## could use some other distribution?
 out <- FocusCH(data,
                fun.cost = lr_Focus, # MdFOCuS for Gaussian Model(change in mean)
                common_difference_step = 1,
                common_ratio_step = 2,
                first_step_qhull = ncol(data) + 50)$opt.cost
 -min(cummin(out))  ### get max stat
 ### pos <- which(diff(res) !=0)+1 ## Store possible stopping times
 ### list(pos=pos, val=res[pos], max_n=-min(out[1:n]))
}

#------------------------------------------------------------------------|
#' @title sim_Max_ChangeInMeanAndVar
#'
#' @description Function to simulate 1d data and to return the max statistics of MdFOCuS for Gaussian Model(change in mean and variance)
#' @param n number of data points
#' @param sim.data distribution (by default, rnorm)
#' @return max statistics
#'
sim_Max_ChangeInMeanAndVar <- function(n,
                                       sim.data = rnorm,
                                       min.size = 2,
                                       minV = 0.001,
                                       col_index = 1:3) {
 x <- sim.data(n)
 data <- cbind(x, x^2) ## could use some other distribution?
 out <- FocusCH_var(data,
                    fun.cost = lr_Focus_var,  # MdFOCuS for Gaussian Model(change in mean and variance)
                    min.size = min.size,
                    minV = minV,
                    common_difference_step = 1,
                    common_ratio_step = 2,
                    first_step_qhull = ncol(data) + 50,
                    col_index = col_index)$opt.cost
 -min(cummin(out))  ### get max stat
}

################################################################################
########################### END ################################################
################################################################################

