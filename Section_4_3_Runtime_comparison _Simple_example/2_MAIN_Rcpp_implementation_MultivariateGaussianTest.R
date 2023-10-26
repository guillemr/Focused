################################################################################
##    SECTION 4.3.1 Runtime of MdFOCuS by R/C++ package "focus"               ##
##    p-variate Independent Gaussian Model without change (p =2,3)            ##
################################################################################
# Test description--------------------------------------------------------|
# Using RR/C++ package "focus" (based on the qhull library)
# For change in mean of a p-variate Gaussian case R/C++ package "focus".  |      
#-------------------------------------------------------------------------|

# Remark------------------------------------------------------------------|                             
# on an Intel(R) Xeon(R) Gold 6252 CPU @ 2.10GHz processor                |                                  
#-------------------------------------------------------------------------|


rm(list = ls())
#packages----------------------------------------------------------------|
#install.packages("remotes")
library(remotes)
remotes::install_github("lpishchagina/focus", force = TRUE)
library(focus)
#parameters---------------------------------------------------------------------
N_<- 5*10^6
P_ <- 3
Cost_ <- "gauss"
mean <- 0
Method_ <- "FOCuS"
Alpha_ <- 2
Beta_ <- 1

set.seed(21)

ts_ <- generate_ts(type = Cost_, p = P_, n = N_, changes = NULL, means = matrix(mean, ncol = 1, nrow = P_))

# Recover runtime for 1 to p dimensions-----------------------------------------
sapply(1:P_, FUN=function(i){
  system.time(res_Mean <- getChangePoints(ts_[, 1:i, drop = F], 
                                          method = Method_,
                                          cost = Cost_,
                                          common_difference_step = Beta_,
                                          common_ratio_step = Alpha_)
  )[3]
}
)

## on an Intel(R) Xeon(R) Gold 6252 CPU @ 2.10GHz processor
## elapsed elapsed elapsed 
##  p=1: 48.363 < 1min  p=2: 238.677 <4min p=3 :1139.583 <19min
