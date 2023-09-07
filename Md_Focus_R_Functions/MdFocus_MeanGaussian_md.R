################################################################################
##           MdFOCuS Implementation by R package "geometry"                   ##
##     Multivariate Independent Gaussian Model(change in mean)                ##
################################################################################

# Remark------------------------------------------------------------------|
#
#-------------------------------------------------------------------------|


#packages----------------------------------------------------------------|
library(geometry)

#functions---------------------------------------------------------------|
#------------------------------------------------------------------------|
#' @title cost_Focus0
#'
#' @description Function to calculate of functional cost(the l2-cost) of MdFOCuS0 (pre-change mean is known and is equal 0)
#'  given the left_cumsum matrix of all considered changepoint candidates
#' @param left_cumsum cumsum of data at considering time
#' @param sum_squares cnst for cost at considering time
#' @return cost of all changepoint candidates

cost_Focus0 <- function(left_cumsum, sum_squares) {
  cand <- 1:(nrow(left_cumsum) - 1)
  right_cumsum <- scale(left_cumsum, left_cumsum[nrow(left_cumsum), ], scale = F)
  -rowSums(right_cumsum[cand, 2:ncol(right_cumsum), drop = FALSE]^2) / abs(right_cumsum[cand, 1]) + sum_squares
}

#------------------------------------------------------------------------|
#' @title lr_Focus0
#'
#' @description Function to return the l2-cost of no change (assuming first segment mean = 0) minus the cost of change given the left_cumsum matrix of all considered changes (- the likelihood ratio test)
#' @param left_cumsum leftcumsum of
#' @param sum_squares
#' @return cost of all considered changes

lr_Focus0 <- function(left_cumsum, sum_squares) {
  cand <- 1:(nrow(left_cumsum) - 1)
  right_cumsum <- scale(left_cumsum, left_cumsum[nrow(left_cumsum), ], scale = F)
  -rowSums(right_cumsum[cand, 2:ncol(right_cumsum), drop = FALSE]^2) / abs(right_cumsum[cand, 1])
}

#------------------------------------------------------------------------|
#' @title cost_Focus
#'
#' @description Function to return the l2-cost given the left_cumsum matrix of all considered changes
#' @param left_cumsum
#' @param sum_squares
#' @return cost of all considered changes

cost_Focus <- function(left_cumsum, sum_squares){
  cand <- 1:(nrow(left_cumsum) - 1)
  right_cumsum <- scale(left_cumsum, left_cumsum[nrow(left_cumsum), ], scale=F)
  ncand <- nrow(left_cumsum) - 1
  - rowSums(right_cumsum[cand, 2:ncol(right_cumsum), drop = FALSE]^2) / abs(right_cumsum[cand, 1]) -
  rowSums(left_cumsum [cand, 2:ncol(left_cumsum), drop = FALSE]^2) / abs(left_cumsum[cand, 1]) + sum_squares
}

#------------------------------------------------------------------------|
#' @title lr_Focus
#'
#' @description Function to return the l2-cost of no change  minus the cost of change given the left_cumsum matrix of all considered changes (- the likelihood ratio test)
#' @param left_cumsum
#' @param sum_squares
#' @return cost of all considered changes

lr_Focus <- function(left_cumsum, sum_squares){
  cand <- 1:(nrow(left_cumsum)-1)
  right_cumsum <- scale(left_cumsum, left_cumsum[nrow(left_cumsum), ], scale=F)
  ncand <- nrow(left_cumsum)-1
  -rowSums(right_cumsum[cand, 2:ncol(right_cumsum), drop=FALSE]^2) / abs(right_cumsum[cand, 1]) -
  rowSums(left_cumsum [cand, 2:ncol(left_cumsum) , drop=FALSE]^2) / abs(left_cumsum[cand, 1]) +
  sum(left_cumsum [nrow(left_cumsum), 2:ncol(left_cumsum) , drop=FALSE]^2) / abs(left_cumsum[nrow(left_cumsum), 1])
}

#------------------------------------------------------------------------|
#' @title FocusCH
#'
#' @description  MdFOCuS algorithm for Gaussian model (change in mean)(using the convex hull algorithm of Quickhull)
#' @param data p-variate time series
#' @param fun.cost type of functional cost : cost_Focus0 or cost_Focus.
#' By defult, fun.cost = cost_Focus0 (pre-change mean is known and is equal 0).
#' @param common_difference_step beta parameter (by default, beta = 1)
#' @param common_ratio_step alpha parameter (by default, alpha = 2)
#' @param first_step_qhull the first step for pruning (using convex hull)
#' By default, first_step_qhull = p + 5
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

FocusCH <- function(data,
                    fun.cost = cost_Focus0,
                    common_difference_step = 1,
                    common_ratio_step = 1,
                    first_step_qhull = ncol(data) + 5,
                    threshold = Inf) {
  # Initialization-------------
  ##get data parameters
  n <- nrow(data)
  p <- ncol(data)

  ## get cnsts for costs
  data_cumsum_square <- cumsum(rowSums(data^2))

  ##get points P(i), i in {1,n}
  data <- cbind (1, data)
  data_left_cumsum  <- apply(data, 2, cumsum)

  ##pruning step
  next_step_qhull <- first_step_qhull

  ## Output candidate list
  list_cand <- list(index = integer(n), # set of candidates
                  nb = 0, #number of candidates
                  nb_at_step = integer(n), #vector of the number of candidates at each iteration
				          opt.cost = numeric(n),   #vector of the functional cost at each iteration
				          opt.change = integer(n))  #vector of the optimal changepoint candidate at each iteration

  ##First iteration
  list_cand$nb                  <- list_cand$nb + 1
  list_cand$index[list_cand$nb] <- 1
  list_cand$nb_at_step[1]       <- 1
  #------------------------------------
  #For loop----------------------------
  for (i in 2:n) {  # start at 2
    ## add a candidate
    list_cand$nb                  <- list_cand$nb + 1;
    list_cand$index[list_cand$nb] <- i;
    list_cand$nb_at_step[i]       <- list_cand$nb

    ##First step: search of optimal changepoint candidate (index and maximization)
    index_cand <- list_cand$index[1:list_cand$nb]
    left_mean  <- data_left_cumsum[index_cand, ] #get cumsum for cost at time i
    cost_cand               <- fun.cost(left_mean, data_cumsum_square[i]) #get all costs at time i
    list_cand$opt.cost[i]   <- min(cost_cand) # get optimal cost
    list_cand$opt.change[i] <- index_cand[which.min(cost_cand)]  # get optimal changepoint candidate

    ## Second step: Pruning by the Quickhull algorithm from the package "geometry"
    if(list_cand$nb >= next_step_qhull | i == n) {
      hull <- convhulln(left_mean) # get convex hull of {P(i)} i in set of candidates
      on_the_hull <- sort(unique(as.vector(hull))) # get vertices
      list_cand$nb <- length(on_the_hull) #get number of vertices
      list_cand$index[1:list_cand$nb] <- index_cand[on_the_hull]; #update the set of candidates
      ## update next_step_qhull
      next_step_qhull <- common_ratio_step * list_cand$nb + common_difference_step
    }

    ## if the (minus) optimal cost is over the threshold end the algorithm
    if (- list_cand$opt.cost[i] >= threshold) {
      break
    }

  }
  return(list_cand)
}

################################################################################
########################### END ################################################
################################################################################

