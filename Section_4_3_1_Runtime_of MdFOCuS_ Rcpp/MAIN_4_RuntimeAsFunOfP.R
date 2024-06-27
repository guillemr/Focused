################################################################################
##                      Empirical Runtime of MdFOCuS by Rcpp                  ##
##    p-variate Gaussian and Poisson Models without change (p =1,..,5)        ##
################################################################################


# Test description--------------------------------------------------------|
# We generate ts ~ N(O,Ip) or P(1p) in dimension p = 1,.., 5              | 
# with n = 2^(10:23) data points (without changes)  and estimate          |
# the runtimes of MdFOCuS and MdFOCuS0 by R/C++ package "focus".          |                                            |    
#-------------------------------------------------------------------------|

# Remark------------------------------------------------------------------|
# The test was performed using parallel computing on the server           |
# The number of cores: mc.cores = 8, where                                |
# nbSimus_ = 100 is the number of simulations.                            |
# The test results are saved in files of the following type:              |
# 'Runtime_alpha_2_beta_1_P_2_FOCuS_gauss_it_100.txt'                     |
# 'Runtime_alpha_2_beta_1_P_2_FOCuS_poisson_it_100.txt'                   |
# Time limit: TimeLimit = 15mn                                            |
#-------------------------------------------------------------------------|

#packages----------------------------------------------------------------|
#install.packages("remotes")
library(remotes)
#remotes::install_github("lpishchagina/focus", force = TRUE)
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
#' @example 
#' Runtime(cnst_=1,N_=1000,P_=2, Method_ = "FOCuS0", Cost_ = 'poisson')
Runtime <- function(cnst_,
                    N_,
                    P_,
                    Method_,
                    Cost_ ,
                    Beta_ = 1,
                    Alpha_ = 2
                    
) {
  set.seed(cnst_ + 10)
  if (Cost_ == 'gauss')  {
    mean <- 0
    first_step_qhull_ <- 1
  }
  else if (Cost_ == "poisson")  {
    mean <- 1
    first_step_qhull_ <- 14
  }
  ts_ <- generate_ts(type = Cost_, p = P_, n = N_, changes = NULL, means = matrix(mean, ncol = 1, nrow = P_))
  res <- system.time(getChangePoints(ts_,
                                     method = Method_,
                                     cost = Cost_,
                                     common_difference_step = Beta_,
                                     common_ratio_step = Alpha_, 
                                     first_step_qhull = first_step_qhull_))[[1]]
  return (res);
}
#Runtime (1,100,2,Method_="FOCuS0",Cost_  = 'poisson', Beta_ = 1, Alpha_ = 2,first_step_qhull_ = 14) 
#method parameters--------------------------------------------------------|
method = c('FOCuS0', 'FOCuS')
cost = c("gauss", "poisson")
#parameters---------------------------------------------------------------|
degree <- 10:23
n_ <-2^(degree)
dim_ <- c(1, 2, 3, 4, 5)
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
                                          Beta_ = 1,
                                          Alpha_ = 2
          ), mc.cores = 100)) #CHANGE FOR YOUR PC
        Time <- max(TableOfRuntime[index , ])
        index <- index + 1
      }
      #save result : average runtime
      write.table(data.frame (n = n_, rowMeans(TableOfRuntime, na.rm = TRUE)), 
                  paste('Runtime_alpha_2_beta_1_P', dim_[p],method[t],cost[cst],'it',nbSimus_,'.txt',sep = '_'), 
                  row.names = TRUE, col.names = FALSE)
    }
  }
}



################################################################################
########################### END ################################################
################################################################################