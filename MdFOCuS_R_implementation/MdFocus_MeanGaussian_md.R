################################################################################
##           MdFOCuS Implementation by R package "geometry"                   ##
##     Multivariate Independent Gaussian Model(change in mean)                ##
################################################################################

# Remark------------------------------------------------------------------|
#
#-------------------------------------------------------------------------|


#packages----------------------------------------------------------------|
library(geometry)
library(purrr)
library(Rcpp)
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
  # this is just a telescopic sum cancellation trick
  # subtracting from all candidates, the last column, e.g. the total cumulative sum
  right_cumsum <- scale(left_cumsum, left_cumsum[nrow(left_cumsum), ], scale=F)
  ncand <- nrow(left_cumsum)-1
  -rowSums(right_cumsum[cand, 2:ncol(right_cumsum), drop=FALSE]^2) / abs(right_cumsum[cand, 1]) -
  rowSums(left_cumsum [cand, 2:ncol(left_cumsum) , drop=FALSE]^2) / abs(left_cumsum[cand, 1]) +
  sum(left_cumsum [nrow(left_cumsum), 2:ncol(left_cumsum) , drop=FALSE]^2) / abs(left_cumsum[nrow(left_cumsum), 1])
}

# -------------------------------------# 
# get global optimum

get_glo_opt <- function(left_cusum, sum_squares, cost=cost_Focus0) {
  out_costs <- cost(left_cusum, sum_squares)
 out <- list(
   opt.change = which.min(out_costs),
   opt.cost = min(out_costs)
 )
 out
}

#------------------------------------------------------------------------|
#' cost_lr_partial
#'
#' @description Function for computing the likelihood ratio test for increasing density regimes, e.g. change in 1-d, 2-d, ...,
#' 
#' @param left_cumsum 
#' @param sum_squares 
#' @param null_known 
#'
#' @return
#' @export
#'
#' @examples
cost_lr_partial <- function(left_cumsum, sum_squares, null_known = FALSE){
  
  cand <- 1:(nrow(left_cumsum)-1)
  right_cumsum <- scale(left_cumsum, left_cumsum[nrow(left_cumsum), ], scale=F)
  ncand <- nrow(left_cumsum)-1
  right_cost <- -(right_cumsum[cand, 2:ncol(right_cumsum), drop=FALSE]^2) / abs(right_cumsum[cand, 1]) 
  
  if(null_known) {
    indep_cusum <-  right_cost
  } else {
    left_cost <- - (left_cumsum [cand, 2:ncol(left_cumsum) , drop=FALSE]^2) / abs(left_cumsum[cand, 1]) 
    null_cost <-  ((left_cumsum [nrow(left_cumsum), 2:ncol(left_cumsum) , drop=FALSE]^2) / abs(left_cumsum[nrow(left_cumsum), 1]))[rep(1, ncand), ] # possible bug
    indep_cusum <- right_cost + left_cost + null_cost
  }
  
  #apply(indep_cusum, 1, sort) |> apply(2, cumsum) |> t()
  indep_cusum
}

cost_lr_partial0 <- \(...) cost_lr_partial(..., null_known = T)


## the min and argminum
cppFunction('
List sortCumsumSelectedMinArgmin_(NumericMatrix mat, IntegerVector columns) {
    int nrow = mat.nrow();
    int ncol = mat.ncol();
    int ncol_output = columns.size();
    NumericVector min_vals(ncol_output);
    IntegerVector argmin_vals(ncol_output);
    
    for (int i = 0; i < nrow; i++) {
        std::vector<double> row(ncol);
        
        // Copy the row to a std::vector
        for (int j = 0; j < ncol; j++) {
            row[j] = mat(i, j);
        }
        
        // Sort the std::vector
        std::sort(row.begin(), row.end());
        
        // Calculate cumulative sum
        for (int j = 1; j < ncol; j++) {
            row[j] += row[j - 1];
        }
        
        // Copy the specified columns of the cumulative sum back to the result matrix
        for (int k = 0; k < ncol_output; k++) {
            if(row[columns[k] - 1] < min_vals[k]){
            min_vals[k] = row[columns[k] - 1];
            argmin_vals[k] = i;
            }
        }
    }
    
    
    
    return List::create(
        Named("min_vals") = min_vals,
        Named("argmin_vals") = argmin_vals
    );
}
')

#------------------------------------------------------------------------|
#' get_partial_opt
#'
#' @description
#' Get the optimal values across increasing density regimes
#' 
#' @param left_cusum 
#' @param sum_squares 
#' @param cost A cost function, like cost_lr_partial, should output a matrix of dimensions #T x p, where T is the set of opt changepoints considered
#'
#' @return
#' @export
#'
#' @examples
get_partial_opt <- function(left_cusum, sum_squares, cost=cost_lr_partial, which_par = ncol(left_cusum) - 1) {
  out_costs <- cost(left_cusum, sum_squares)

  sum_of_the_max <- out_costs |> apply(2, min) |> sum()

  #exact_partials <- apply(out_costs, 1, sort) |> apply(2, cumsum) |> t()

  out_cpp <- sortCumsumSelectedMinArgmin_(out_costs, c(1, which_par))

  out <- list(
    opt.change = out_cpp$argmin_vals,
    opt.cost = c(m = out_cpp$min_vals[1],
                 sm = sum_of_the_max,
                 ex = out_cpp$min_vals[2:length(out_cpp$min_vals)])
  )

  # out <- list(
  #   opt.change = exact_partials[, which_par, drop=F] %>% apply(2, which.min),
  #   opt.cost = c(m = min(exact_partials[ , 1]),
  #                sm = sum_of_the_max,
  #                ex = exact_partials[, which_par, drop=F] %>% apply(2, min))
  # )

  # out <- list(
  #   opt.change = apply(out_costs, 2, which.min),
  #   opt.cost = apply(out_costs, 2, min)
  # )
  out
}

get_3_opts_reg <- function(left_cusum, sum_squares, cost=cost_lr_partial) {
  out_costs <- cost(left_cusum, sum_squares)

  partial_sums <- out_costs |> apply(2, min) |> sort() |> cumsum()
  exact <- rowSums(out_costs)

  out <- list(
    opt.change = exact |> which.min(),
    opt.cost = c(m = partial_sums[1],
                 sm = partial_sums[length(partial_sums)],
                 ex = exact |> min())
  )
  out
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
                    get_opt_cost = get_glo_opt,
                    common_difference_step = 1,
                    common_ratio_step = 2,
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
                    opt.cost = vector(mode = "list", length = n),   #vector of the functional cost at each iteration
                    opt.change = vector(mode = "list", length = n))  #vector of the optimal changepoint candidate at each iteration
  
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
    
    out <- get_opt_cost(left_mean, data_cumsum_square[i]) # get optimal cost
    list_cand$opt.change[[i]] <- out$opt.change
    list_cand$opt.cost[[i]] <- out$opt.cost
    
    
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
    if(length(list_cand$opt.cost[[i]]) == length(threshold))
      change_detected <- sum(- list_cand$opt.cost[[i]] >= threshold)
    else
      warning("The number of statistics to test does not match the length of the threshold")
    if (change_detected) {
      break
    }
    
  }
  return(list_cand)
}


FocusCH_HighDim <- function(data,
                    get_opt_cost = get_glo_opt,
                    common_difference_step = 1,
                    common_ratio_step = 2,
                    first_step_qhull = ncol(data) + 5,
                    threshold = Inf,
                    dim_indexes = map2(seq(0, ncol(data)-1, by=2), seq(1, ncol(data), by=2), \(f, s) c(f, s %% ncol(data)) + 1)) {
  # Initialization-------------
  ##get data parameters
  n <- nrow(data)
  p <- ncol(data)
  
  ## get cnsts for costs
  data_cumsum_square <- cumsum(rowSums(data^2))
  
  ##get points P(i), i in {1,n}
  data_left_cumsum  <- apply(cbind (1, data), 2, cumsum)
  
  ##pruning step
  next_step_qhull <- first_step_qhull
  
  ## Output candidate list
  list_cand <- list(index = integer(n), # set of candidates
                    nb = 0, #number of candidates
                    nb_at_step = integer(n), #vector of the number of candidates at each iteration
                    opt.cost = vector(mode = "list", length = n),   #vector of the functional cost at each iteration
                    opt.change = vector(mode = "list", length = n))  #vector of the optimal changepoint candidate at each iteration
  
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
    
    out <- get_opt_cost(left_mean, data_cumsum_square[i]) # get optimal cost
    list_cand$opt.change[[i]] <- out$opt.change
    list_cand$opt.cost[[i]] <- out$opt.cost
    
    
    ## Second step: Pruning by the Quickhull algorithm from the package "geometry"
    if(list_cand$nb >= next_step_qhull | i == n) {
      on_the_hull <-
        map(dim_indexes, \(i) convhulln(left_mean[, c(1, i + 1)]) |>
                  as.vector() |>
                  unique()) |>
        unlist() |>
        unique() |>
        sort()
      
      
      #hull <- convhulln(left_mean) # get convex hull of {P(i)} i in set of candidates
      #on_the_hull <- sort(unique(as.vector(hull))) # get vertices
      list_cand$nb <- length(on_the_hull) #get number of vertices
      list_cand$index[1:list_cand$nb] <- index_cand[on_the_hull]; #update the set of candidates
      ## update next_step_qhull
      next_step_qhull <- common_ratio_step * list_cand$nb + common_difference_step
    }
    
    ## if the (minus) optimal cost is over the threshold end the algorithm
    if(length(list_cand$opt.cost[[i]]) == length(threshold))
      change_detected <- sum(- list_cand$opt.cost[[i]] >= threshold)
    else
      warning("The number of statistics to test does not match the length of the threshold")
    if (change_detected) {
      break
    }
    
  }
  return(list_cand)
}

FocusCH_HighDim_OPT <- function(data,
                    get_opt_cost = get_glo_opt,
                    common_difference_step = 1,
                    common_ratio_step = 2,
                    first_step_qhull = ncol(data) + 5,
                    threshold = Inf,
                    dim_indexes = map2(seq(0, ncol(data)-1, by=2), seq(1, ncol(data), by=2), \(f, s) c(f, s %% ncol(data)) + 1)) {
  
  update_hull <- function(l_c) {
    if (l_c$nb >= l_c$next_step_qhull | i == n) {
      index_cand <- l_c$index[1:l_c$nb]
      local_left_mean <- data_left_cumsum[index_cand, c(1, l_c$i_d)]
      hull <- convhulln(local_left_mean) # get convex hull of {P(i)} i in set of candidates
      on_the_hull <- sort(unique(as.vector(hull))) # get vertices
      l_c$nb <- length(on_the_hull) # get number of vertices
      l_c$index[1:l_c$nb] <- index_cand[on_the_hull] # update the set of candidates
      ## update next_step_qhull
      l_c$next_step_qhull <- common_ratio_step * l_c$nb + common_difference_step
    }
    return(l_c)
  }

  # Initialization-------------
  ##get data parameters
  n <- nrow(data)
  p <- ncol(data)
  
  ## get cnsts for costs
  data_cumsum_square <- cumsum(rowSums(data^2))
  
  ##get points P(i), i in {1,n}
  data_left_cumsum  <- apply(cbind (1, data), 2, cumsum)
  

  list_cand <- map(dim_indexes, \(i_d) list(index = integer(n), # set of candidates
                                            nb = 0, #number of candidates
                                            nb_at_step = integer(n), #vector of the number of candidates at each iteration
                                            i_d = i_d + 1,
                                            next_step_qhull = first_step_qhull)
  )
  
  opt_list <- list(
    opt.cost = vector(mode = "list", length = n),   #vector of the functional cost at each iteration
    opt.change = vector(mode = "list", length = n)  #vector of the optimal changepoint candidate at each iteration
    )


  ##First iteration
  list_cand <- map(list_cand, \(l_c) {
    l_c$nb                  <- l_c$nb + 1
    l_c$index[l_c$nb] <- 1
    l_c$nb_at_step[1]       <- 1
    l_c
  })

  
  #------------------------------------
  #For loop----------------------------
  for (i in 2:n) {  # start at 2
    ## add a candidate
    list_cand <- map(list_cand, \(l_c) {
      l_c$nb                  <- l_c$nb + 1
      l_c$index[l_c$nb]       <- i
      l_c$nb_at_step[i]       <- l_c$nb
      l_c
    })
    
    ##First step: search of optimal changepoint candidate (index and maximization)

    glo_index_cand <- map(list_cand, \(l_c) l_c$index[1:l_c$nb]) |> 
                        unlist() |>
                        unique()
 
    glo_left_mean  <- data_left_cumsum[glo_index_cand, ] #get cumsum for cost at time i
    out <- get_opt_cost(glo_left_mean, data_cumsum_square[i]) # get optimal cost
    opt_list$opt.change[[i]] <- out$opt.change
    opt_list$opt.cost[[i]] <- out$opt.cost

    # for debug reasons only
    # opt_list$opt.change[[i]] <- rep(0, length(threshold))
    # opt_list$opt.cost[[i]] <- rep(0, length(threshold))

    ## Second step: Pruning by the Quickhull algorithm from the package "geometry"
    list_cand <- map(list_cand, update_hull)

    
    ## if the (minus) optimal cost is over the threshold end the algorithm
    if(length(opt_list$opt.cost[[i]]) == length(threshold)) {
      change_detected <- sum(- opt_list$opt.cost[[i]] >= threshold)
    } else {
      warning("The number of statistics to test does not match the length of the threshold")
    }

    if (change_detected) {
      break
    }
    
  }
  return(opt_list)
}

################################################################################
########################### END ################################################
################################################################################

if(F) {      
  source("Section_4_3_2_OCDlike_Simulation/helper_functions.R")
    
  y_nc <- generate_sequence(n = 1000, p = 3, cp = 999, magnitude = 1, dens = 1, seed = 123) |> t()

  # with the unkown loss 
  system.time(out <- FocusCH(y_nc, get_opt_cost = \(...) get_glo_opt(..., cost=lr_Focus0), threshold = Inf))
  - head(out$opt.cost |> unlist())
}

# testing stuff
if (F) {
  library(tidyverse)
  library(future)
  library(furrr)
  plan(multisession, workers=18)
  
  source("Section_4_3_2_OCDlike_Simulation/helper_functions.R")
  
  y_nc <- generate_sequence(n = 5000, p = 3, cp = 4999, magnitude = 0, dens = 0, seed = 42) |> t()

  # with the known loss  
  cost_lr_partial0 <- \(...) cost_lr_partial(..., null_known = T)
  
  out <- FocusCH(y_nc, get_opt_cost = \(...) get_partial_opt(..., cost = cost_lr_partial0), threshold = rep(Inf, 3))
  out$opt.cost

  # with the unkown loss 
  out <- FocusCH(y_nc, get_opt_cost = get_partial_opt, threshold = rep(Inf, 3))
  out$opt.cost
  
  # tuning the thresholds
  get_one_run <- function(seed, n = 100, p = 3, v_density = 1:p) {
    y_nc <- generate_sequence(n = n, p = p, cp = n-1, magnitude = 0, dens = 0, seed = seed) |> t()
    out <- FocusCH(y_nc, get_opt_cost = get_partial_opt, threshold = rep(Inf, p))
    -reduce(out$opt.cost, rbind)[, v_density] |> apply(2, max)
  }
  
  out1 <- future_map(1:100, get_one_run, p = 3) |> reduce(rbind)
  apply(out1, 2, quantile, prob=exp(-1))
  
  
  ## unknown
  get_one_run <- function(seed, n = 100, p = 3, v_density = 1:p) {
    y_nc <- generate_sequence(n = n, p = p, cp = n-1, magnitude = 0, dens = 0, seed = seed) |> t()
    out <- FocusCH(y_nc, get_opt_cost = \(...) get_partial_opt(..., cost = cost_lr_partial0), threshold = rep(Inf, p))
    -reduce(out$opt.cost, rbind)[, v_density] |> apply(2, max)
  }
  
  out2 <- future_map(1:100, get_one_run, p = 3) |> reduce(rbind)
  apply(out2, 2, quantile, prob=exp(-1))
}


# testing high dimentional

if (F) {
  library(tidyverse)
  library(future)
  library(furrr)
  plan(multisession, workers=18)
  
  source("Section_4_3_2_OCDlike_Simulation/helper_functions.R")
  
  p <- 100
  N <- 2000
  y = generate_sequence(n = N, p = p, cp = 500, magnitude = 0, dens = 0, seed = 42)  

  # focus
  res <- FOCuS_multi_JMLR(y, c(Inf, Inf), mu0 = rep(0, p))
  data.frame(max = max(res$maxs), sum = max(res$sums))

  data <- t(y) # trasposing as the current prototype reads nxp (rather then pxn)
  system.time(res <- FocusCH_HighDim(data, get_opt_cost = \(...) get_partial_opt(..., cost=cost_lr_partial0, which_par = c(5, 25, 100)), common_ratio_step = 2, threshold = rep(Inf, 5)))
  - (res$opt.cost |> reduce(rbind)) |> apply(2, max)

  ocd_det <- ocd_known(c(Inf, Inf, Inf), rep(0, p), rep(1, p))
  system.time(r <- ocd_detecting(y, ocd_det))
}

# testing fast high dimentional

if (F) {
  library(tidyverse)
  library(future)
  library(furrr)
  library(Rcpp)
  plan(multisession, workers=18)
  
  source("exemple_Rcpp.R")

  source("Section_4_3_2_OCDlike_Simulation/helper_functions.R")
  
  p <- 50
  N <- 5000
  y = generate_sequence(n = N, p = p, cp = N-1, magnitude = 0, dens = 0, seed = 123)



  data <- t(y) # trasposing as the current prototype reads nxp (rather then pxn)
  #res <- FocusCH_HighDim_OPT(data, get_opt_cost = \(...) get_partial_opt(..., cost=cost_lr_partial0, which_par = c(5, 25, 100)), common_ratio_step = 2, threshold = rep(Inf, 5))
  system.time(res1 <- FocusCH_HighDim_OPT(data, get_opt_cost = \(...) get_partial_opt(..., cost=cost_lr_partial0, which_par = c(5, 25, 100)), common_ratio_step = 2, threshold = rep(Inf, 5)))
  - (res1$opt.cost |> reduce(rbind)) |> apply(2, max)
  system.time(res <- FocusCH_HighDim(data, get_opt_cost = \(...) get_partial_opt(..., cost=cost_lr_partial0, which_par = c(5, 25, 100)), common_ratio_step = 1.5, threshold = rep(Inf, 5)))
  - (res$opt.cost |> reduce(rbind)) |> apply(2, max)

  ocd_det <- ocd_known(c(Inf, Inf, Inf), rep(0, p), rep(1, p))
  system.time(r <- ocd_detecting(y, ocd_det))


  # res3 <- FocusCH_HighDim(data, get_opt_cost = \(...) get_partial_opt(..., cost=cost_lr_partial0, which_par = c(5, 25, 100)),
  #                                         common_ratio_step = 1.2,
  #                                         dim_indexes = as.list(1:ncol(data)),
  #                                         threshold = rep(Inf, 5))
  # - (res3$opt.cost |> reduce(rbind)) |> apply(2, max)


}