################################################################################
##    Empirical Runtime of MdFOCuS0 by Rcpp as the function of alpha          ##
##    p-variate Independent Gaussian Model without change (p =2,3)            ##
################################################################################

# Test description--------------------------------------------------------|
# We generate ts ~ N(O,Ip)  with n = 10^5 data points (without change),   |
# p=2,3 and alpha= c(1+2^(-4:0), 3, 2^(2:7))                              |
# and estimate the runtimes of MdFOCuS0 by R/C++ package "focus".         |      
#-------------------------------------------------------------------------|

# Remark------------------------------------------------------------------|                             
# The test was performed using parallel computing on the server           |
# The number of cores: mc.cores = nbSimus_, where                         |
# nbSimus_ = 100 is the number of simulations.                            |
# The test results are saved in files :                                   |
# 'Alpha_dependence_2_N_1e+05_FOCuS0_gauss_it_100_.txt' and               |
# 'Alpha_dependence_3_N_1e+05_FOCuS0_gauss_it_100_.txt'                   |                                  
#-------------------------------------------------------------------------|

#packages----------------------------------------------------------------|
#install.packages("remotes")
library(remotes)
remotes::install_github("lpishchagina/focus", force = TRUE)
library(focus)
library(parallel)

#function----------------------------------------------------------------|
#' @title Runtime
#'
#' @description time complexity
#' @param cnst parameter for function set.seed()
#' @param N_ length of time series
#' @param P_ dimension
#' @param Method_ method
#' @param Cost_ type of cost function
#' @param First_step (By default, first step + P_) 
#' @param Beta_  By default, Beta_ = 1
#' @param Alpha_  Value of alpha
#'
#' @return returns time complexity for ts (without changes)
Runtime <- function(cnst_,
                    N_=10^5,
                    P_= 2,
                    Method_ = 'FOCuS0',
                    Cost_  = 'gauss',
                    Beta_ = 1,
                    Alpha_,
                    first_step_qhull_ = 1
) {
  set.seed(cnst_ + 10)
  if (Cost_ == 'gauss')  mean <- 0
  else if (Cost_ == "poisson")  mean <- 1
  ts_ <- generate_ts(type = Cost_, p = P_, n = N_, changes = NULL, means = matrix(mean, ncol = 1, nrow = P_))
  res <- system.time(getChangePoints(ts_,
                                     method = Method_,
                                     cost = Cost_,
                                     common_difference_step = Beta_,
                                     common_ratio_step = Alpha_))[[1]]
  return (res);
}

#method parameters--------------------------------------------------------|
Method = 'FOCuS0'
Cost = 'gauss'
#parameters---------------------------------------------------------------|
N <-10^5
P <-  c(2, 3)
Alpha <- c(1+2^(-4:0), 3, 2^(2:7))
#number of simulations
nbSimus_ <- 100

#calculations-------------------------------------------------------------|
for (p in 1 : length(P)) {
  TableOfRuntime <- matrix(NA, length(Alpha),  nbSimus_)
  for (alpha in 1 : length(Alpha)) {
    alpha_ <- Alpha[alpha]
    TableOfRuntime[alpha, ]  <- do.call(cbind, parallel::mclapply(
      1:nbSimus_, function(i) Runtime (cnst_ = i,
                                       N_ = N,
                                       P_= P[p],
                                       Method_ = Method,
                                       Cost_  = Cost,
                                       Beta_ = 1,
                                       Alpha_ = alpha_,
                                       first_step_qhull_ = 1
      ), mc.cores = 100))##CHANGE FOR YOUR PC
  }
  #save result : average run time
  write.table(data.frame (alpha = Alpha, rowMeans(TableOfRuntime, na.rm = TRUE)), 
              paste('Alpha_dependence',P[p],'N', N,Method,Cost,'it',nbSimus_,'.txt',sep = '_'), 
              row.names = TRUE, col.names = FALSE)
}

################################################################################
########################### END ################################################
################################################################################
