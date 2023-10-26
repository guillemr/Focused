################################################################################
##    Empirical Runtime of dyadic MdFOCuS0 by Rcpp as the function of q_min   ##
##    p-variate Independent Gaussian Model without change (p = 1,2,3)         ##
################################################################################

# Test description--------------------------------------------------------|
# We generate ts ~ N(O,Ip)  with n = 10^5 data points (without change),   |
# p= 1, 2,3 and qmin = 8                                                  |
# and estimate the runtimes of dyadic MdFOCuS0
# by R/C++ package "dyadicMdFOCuS".                                       |      
#-------------------------------------------------------------------------|

# Remark------------------------------------------------------------------|                             
# The test was performed using parallel computing on the server           |
# The number of cores: mc.cores = nbSimus_, where                         |
# nbSimus_ = 100 is the number of simulations.                            |
# The test results are saved in files :                                   |
# 'Qmin_dependence_2_N_1e+05_FOCuS0_gauss_it_100_.txt' and               |
# 'Qmin_dependence_3_N_1e+05_FOCuS0_gauss_it_100_.txt'                   |                                  
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
                    N_=10^5,
                    P_= 2,
                    Method_ = 'FOCuS0',
                    Cost_  = 'gauss',
                    q_min = 3
) {
  set.seed(cnst_ + 10)
  mean <- 0
#  if (Cost_ == 'gauss')  mean <- 0
#  else if (Cost_ == "poisson")  mean <- 1
  ts_ <- generate_ts(type = Cost_, p = P_, n = N_, changes = NULL, means = matrix(mean, ncol = 1, nrow = P_))
  res <- system.time(dyadMdFOCuS(ts_, method ='FOCuS0', cost  = 'gauss', qmin = q_min))[[1]]
  return (res);
}


#Runtime(cnst_=21, N_=10^4, P_= 2,Method_ = 'FOCuS0',Cost_  = 'gauss', q_min = 2)

#Runtime(cnst_=21, N_=10^4, P_= 3, Method_ = 'FOCuS0',Cost_  = 'gauss', q_min = 2)

#method parameters--------------------------------------------------------|
Method = 'FOCuS0'
Cost = 'gauss'
#parameters---------------------------------------------------------------|
N <-10^5
P <-  c(1, 2, 3)
#number of simulations
nbSimus_ <- 100


Qmin <- as.integer(log2(P+2))
Qmax <- as.integer(log2(N))

vect_qmin <- 3:12

#calculations-------------------------------------------------------------|
for (p in 1 : length(P)) {
  TableOfRuntime <- matrix(NA, length(vect_qmin),  nbSimus_)
  for (qmin in 1 : length(vect_qmin)) {
    qmin_ <- vect_qmin[qmin]
    TableOfRuntime[qmin, ]  <- do.call(cbind, parallel::mclapply(
      1:nbSimus_, function(i) Runtime (cnst_ = i,
                                       N_ = N,
                                       P_= P[p],
                                       Method_ = Method,
                                       Cost_  = Cost,
                                       q_min =  qmin_ 
      ), mc.cores = nbSimus_))##CHANGE FOR YOUR PC
  }
  #save result : average run time
  write.table(data.frame (qmin = vect_qmin, rowMeans(TableOfRuntime, na.rm = TRUE)), 
              paste('Qmin_dependence',P[p],'N', N,Method,Cost,'it',nbSimus_,'.txt',sep = '_'), 
              row.names = TRUE, col.names = FALSE)
}

################################################################################
########################### END ################################################
################################################################################







