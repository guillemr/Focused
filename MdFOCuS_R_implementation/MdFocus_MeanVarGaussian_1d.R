################################################################################
##           MdFOCuS Implementation by R package "geometry"                   ##
##    Multivariate Independent Gaussian Model(change in mean and variance)    ##
################################################################################

# Remark------------------------------------------------------------------|
#
#-------------------------------------------------------------------------|

#packages----------------------------------------------------------------|

#functions---------------------------------------------------------------|
#------------------------------------------------------------------------|
#' @title getV
#'
#' @description
#' @param mat_cumsum P(i)? cumsum of data at considering time
#' @param minV By default, minV = 1
#' @return

getV <- function(mat_cumsum, minV = 1) {
  data.var <- (abs(mat_cumsum[, 3]) - mat_cumsum[, 2]^2  / abs(mat_cumsum[, 1]) ) / abs(mat_cumsum[, 1])
  est.var <- pmax(data.var, minV)

  abs(mat_cumsum[, 1])*( log(est.var) +  data.var/est.var)#
}

#------------------------------------------------------------------------|
#' @title lr_Focus_var
#'
#' @description Function to return the l2-cost of no change  minus the cost of change given the left_cumsum matrix of all considered changes
#' @param left_cumsum of the considered changes
#' @param all_cumsum cumsum of the whole chunk of data
#' @return cost of all considered changes

lr_Focus_var <- function(left_cumsum, all_cumsum, minV = 1) {
  right_cumsum <- scale(left_cumsum, all_cumsum, scale = F)
  right_v <- getV(right_cumsum[, , drop = FALSE], minV)
  left_v  <- getV(left_cumsum [, , drop = FALSE], minV)
  all_v   <- getV(all_cumsum)
  out <- right_v + left_v - all_v
  #out[is.infinite(out)] <- 0
  out
}

#------------------------------------------------------------------------|
#' @title FocusCH_var
#' @title FocusCH
#'
#' @description  MdFOCuS algorithm for Gaussian model (change in mean and variance) (using the convex hull algorithm of Quickhull)
#' @param data  matrix with 2 columns : p-variate time series x and p-variate time series x^2
#' @param fun.cost type of functional cost : lr_Focus_var (functional cost for  Gaussian model (change in mean and variance)).
#' @param min.size
#' @param minV
#' @param common_difference_step beta parameter (by default, beta = 1)
#' @param common_ratio_step alpha parameter (by default, alpha = 2)
#' @param first_step_qhull the first step for pruning (using convex hull)
#' By default, first_step_qhull = p + 5
#' @param col_index in case only 0 and 1 take 1:2 if not 1:3
#' @return a list of  elements  = (index, nb, nb_at_step, opt.cost, opt.change).
#'
#' \describe{
#' \item{\code{index}}{is the set of candidates.}
#' \item{\code{nb}}{is the number of candidates.}
#' \item{\code{nb_at_step}}{is a vector of the number of candidates at each iteration.}
#' \item{\code{opt.cost}}{is a vector of the functional cost at each iteration.}
#' \item{\code{opt.change}}{is a number of candidates at each iteration (vector).}
#' \item{\code{opt.cost}}{is a cost or -likelihood-ratio vector of the optimal changepoint candidate at each iteration.}
#' }

FocusCH_var <- function(data,
                        fun.cost = lr_Focus_var,
                        min.size = 2,
                        minV = 1,
                        common_difference_step = 1,
                        common_ratio_step = 2,
                        first_step_qhull = ncol(data) + 5,
                        col_index = 1:3) {
  # Initialization-------------
  ##get data parameters
  n <- nrow(data)
  p <- ncol(data)

  ##get points P(i), i in {1,n}
  data     <- cbind(1, data)
  data_left_cumsum  <- apply(data, 2, cumsum)

  ##pruning step
  next_step_qhull <- first_step_qhull

  ## Output candidate list
  list_cand <- list(index = integer(n), # set of candidates
                    nb = 0, #number of candidates
                    nb_at_step = integer(n), #vector of the number of candidates at each iteration
                    opt.cost = numeric(n),   #vector of the functional cost at each iteration
                    opt.change = integer(n))  #vector of the optimal changepoint candidate at each iteration

  list_cand$nb                  <- 0 #@ Luda: why is it not 1?
  list_cand$nb_at_step[1]       <- 0 #@ Luda: why is it not 1?

  #------------------------------------
  #For loop----------------------------
  for(i in (2*min.size):n) {  ## start at 2
    ## add a candidate
    list_cand$nb                  <- list_cand$nb +1;
    list_cand$index[list_cand$nb] <- i - min.size;
    list_cand$nb_at_step[i]       <- list_cand$nb

    ##First step: search of optimal changepoint candidate (index and maximization)
    index_cand <- list_cand$index[1:list_cand$nb]
    left_cumsum  <- data_left_cumsum[index_cand, , drop = FALSE] #get cumsum for cost at time i
    cost_cand               <- fun.cost(left_cumsum, data_left_cumsum[i, , drop = FALSE], minV = minV) #get all costs at time i
    list_cand$opt.cost[i]   <- min(cost_cand) # get optimal cost
    list_cand$opt.change[i] <- index_cand[which.min(cost_cand)] # get optimal changepoint candidate

    ## Second step: Pruning by the Quickhull algorithm from the package "geometry"
    if(list_cand$nb >= next_step_qhull | i == n) {
      hull <- convhulln(left_cumsum[, col_index, drop = F]) # get convex hull of {P(i)} i in set of candidates
      on_the_hull <- sort(unique(as.vector(hull))) # get vertices
      list_cand$nb <- length(on_the_hull) #get number of vertices
      list_cand$index[1:list_cand$nb] <- index_cand[on_the_hull]; #update the set of candidates
      ## update next_step_qhull
      next_step_qhull <- common_ratio_step * list_cand$nb + common_difference_step
    }
  }
  return(list_cand)
}


################################################################################
########################### END ################################################
################################################################################
