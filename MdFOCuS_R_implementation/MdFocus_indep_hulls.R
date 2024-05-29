FocusCH_HighDim_ind <- function(data,
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