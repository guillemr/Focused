################################################################################
##               Empirical Runtime of dyadic MdFOCuS by Rcpp                  ##
##                       p-variate Gaussian without change (p =1,2,3)         ##
################################################################################


# Test description--------------------------------------------------------|
# We generate ts ~ N(O,Ip) or P(1p) in dimension p = 1,2,3                | 
# with n = 2^(10:23) data points (without changes)  and estimate          |
# the runtimes of MdFOCuS and MdFOCuS0 by R/C++ package "focus".          |                                            |    
#-------------------------------------------------------------------------|

# Remark------------------------------------------------------------------|
# The test was performed using parallel computing on the server           |
# The number of cores: mc.cores = 8, where                                |
# nbSimus_ = 100 is the number of simulations.                            |
# The test results are saved in files of the following type:              |
# 'Runtime_alpha_2_beta_1_P_2_FOCuS_gauss_it_100.txt'                     |
# Time limit: TimeLimit = 3mn                                             |
#-------------------------------------------------------------------------|

#packages----------------------------------------------------------------|
#install.packages("remotes")
library(remotes)
remotes::install_github("lpishchagina/dyadicMdFOCuS", force = TRUE)
library(dyadicMdFOCuS)
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
#' @param qmin 
#'
#' @return returns time complexity for ts (without changes)
Runtime <- function(cnst_,
                    N_,
                    P_,
                    Method_ = 'FOCuS0',
                    Cost_  = 'gauss',
                    q_min = 8
) {
  set.seed(cnst_ + 10)
  mean <- 0
  #  if (Cost_ == 'gauss')  mean <- 0
  #  else if (Cost_ == "poisson")  mean <- 1
  ts_ <- generate_ts(type = Cost_, p = P_, n = N_, changes = NULL, means = matrix(mean, ncol = 1, nrow = P_))
  res <- system.time(dyadMdFOCuS(ts_, method ='FOCuS0', cost  = 'gauss', qmin = q_min))[[1]]
  return (res);
}

#method parameters--------------------------------------------------------|
method = c('FOCuS0', 'FOCuS')
cost = c("gauss")
#parameters---------------------------------------------------------------|
degree <- 10:23
n_ <-2^(degree)
dim_ <- c(1, 2, 3)
qmin_opt <- c(6, 7, 8)
#number of simulations
nbSimus_ <- 100
#Time limit
TimeLimit <- 900
#calculations-------------------------------------------------------------|
for (p in 1 : length(dim_)) {
  for ( t in 1 : length(method)) {
    for ( cst in 1 : length(cost)) {
      TableOfRuntime <- matrix(NA, length(n_),  nbSimus_)
      index <- 1
      Time <- 0
      while ((index <= length(n_)) && ((2*Time) <= TimeLimit)) {
        N <- n_[index]
        TableOfRuntime[index , ]  <- do.call(cbind, parallel::mclapply(
          1:nbSimus_, function(i) Runtime(cnst_ =i,
                                          N_ = N,
                                          P_ = dim_[p],
                                          Method_ = method[t],
                                          Cost_ = cost[cst],
                                          q_min = qmin_opt[p]
          ), mc.cores = nbSimus_)) #CHANGE FOR YOUR PC
        Time <- max(TableOfRuntime[index , ])
        index <- index + 1
      }
      #save result : average runtime
      write.table(data.frame (n = n_, rowMeans(TableOfRuntime, na.rm = TRUE)), 
                  paste('Runtime_qmin',qmin_opt[p] ,'P', dim_[p],method[t],cost[cst],'it',nbSimus_,'.txt',sep = '_'), 
                  row.names = TRUE, col.names = FALSE)
    }
  }
}



################################################################################
########################### END ################################################
################################################################################