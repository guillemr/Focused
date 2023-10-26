################################################################################
##         	                       ocd runtime                                ##
##                              dimension p=3                                 ##
################################################################################
#packages----------------------------------------------------------------|
library(ocd)
#source----------------------------------------------------------------|
source("helper_functions_ocd.R")
#parameters---------------------------------------------------------------|
p <- 3
N <- 5*1e5
REPS <- 300
target_arl <- N/2

#data---------------------------------------------------------------------------
# data with no change (for testing the empirical avg run length)
Y_nc <- lapply(1:REPS, function(i) generate_sequence(n = N, p = p, cp = 500, magnitude = 0, dens = 0, seed = i))

# data to estimate the mu0 value
Y_train <- lapply(1:REPS, function(i) generate_sequence(n = 500, p = p, cp = 199, magnitude = 0, dens = 0, seed = 600 + i))

# data to train the monte-carlo treshold
Y_monte_carlo <- lapply(1:REPS, function(i) generate_sequence(n = target_arl + 100, p = p, cp = 500, magnitude = 0, dens = 0, seed = i))


# threshold calculation---------------------------------------------------------
### change detection (no-stop for runtime)
ocd_detecting <- function (Y, detector) {
  for (t in seq_len(ncol(Y))){
    detector <- getData(detector, Y[,t])
    #if (is.numeric(attr(detector, "status")))
    #  break
  }
  list(t = t, det = detector)
}


system.time(ocd_stat <- MC_ocd_v6(Y_monte_carlo, 1, "auto", training_data = Y_train))


q_prob <- .68
ocd_est_thres <- apply(ocd_stat, 2, quantile, prob = q_prob)

ocd_est_thres

#diag    off_d    off_s 
#16.54134 30.09374 30.05722 
#Average run time---------------------------------------------------------------
i <- 1
#data generation
set.seed(21)
y <- t(matrix(rnorm(p*N), nr=N, nc=p))

y_tr <- Y_train[[i]]
ocd_det <- ocd_training(y_tr, ocd_est_thres)

system.time(r <- ocd_detecting(y, ocd_det))
## on an Intel(R) Xeon(R) Gold 6252 CPU @ 2.10GHz processor
## elapsed elapsed elapsed 
#user  system elapsed 
#195.853   2.732 198.583 
################################################################################
########################### END ################################################
################################################################################
