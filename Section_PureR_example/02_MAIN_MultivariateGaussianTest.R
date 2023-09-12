#-----------------------------------------------------------------------------|
# SECTION 2.4 Using Simple R Implementation
# With Geometry R package 
# For change in mean of a multivariate Gaussian                        |

rm(list = ls())
#packages-----------------------------------------------------------------|
library(geometry)

# LOADING Md FOCUS IN R FUNCTIONS
source("../Md_Focus_R_Functions/MdFocus_MeanGaussian_md.R")                # MdFOCuS for Gaussian Model(change in mean)

n <- 5*10^5
p <- 3

set.seed(21)
data1 <- matrix(rnorm(p*n), nr=n, nc=p)

## Recover runtime for 1 to p dimensions
##
sapply(1:p, FUN=function(i){
		system.time(res_Mean <- FocusCH(data1[, 1:i, drop=F], 
                    fun.cost = lr_Focus,
                    common_difference_step = 1, 
                    common_ratio_step = 2, 
                    first_step_qhull = ncol(data1) + 5)
		)[3]
		}
)

## on an Intel(R) Xeon(R) Gold 6252 CPU @ 2.10GHz processor
## elapsed elapsed elapsed 
## 47.259  54.103 145.723 
                    






