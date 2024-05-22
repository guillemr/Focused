##install.packages("remotes", "ocd")
# library(remotes)
##remotes::install_github("gtromano/FOCuS", force = TRUE)
library(ocd)
library(FOCuS)
#library(parallel)
#library(tidyverse)
library(furrr)
library(dplyr)

generate_sequence <- function(n = 1000, p = 100, cp = 200, sd = 1, magnitude = 1, dens = .1, seed = 42){
  
  set.seed(seed)
  noise <- matrix(rnorm(p * n, sd = sd), nr = p, nc = n)
  
  pre_change <- rep(0, p)
  #pre_change <- rnorm(p, mean = 0)
  
  #set.seed(5000 + magnitude * seed)
  zi <- rnorm(floor(p * dens))
  post_means <- magnitude * zi / sqrt(sum(zi^2))
  post_change <- pre_change + c(post_means, rep(0, floor(p * (1-dens))))
  noise + cbind(matrix(pre_change, nr = p, nc = cp), matrix(post_change, nr = p, nc = n - cp))
  
}


# this function takes treshold and known values for pre and post change mean and return an instance of
# an ocd detector ready to do some changepoint monitoring.
ocd_known <- function (thresh, mu0, sd0){
  detector <- ChangepointDetector(dim=length(mu0), method='ocd', beta=1, thresh=thresh)
  detector <- setBaselineMean(detector, mu0)
  detector <- setBaselineSD(detector, sd0)
  
  setStatus(detector, 'monitoring')
}


# this function takes some training observations and return an instance of
# an ocd detector ready to do some changepoint monitoring.
ocd_training <- function (Y_train, thresh){
  detector <- ChangepointDetector(dim=nrow(Y_train), method='ocd', beta=1, thresh=thresh)
  detector <- setStatus(detector, 'estimating')
  for (t in seq_len(ncol(Y_train))) {
    detector <- getData(detector, Y_train[,t])
  }
  setStatus(detector, 'monitoring')
}

ocd_detecting <- function (Y, detector) {
  for (t in seq_len(ncol(Y))){
    detector <- getData(detector, Y[,t])
    if (is.numeric(attr(detector, "status")))
      break
  }
  list(t = t, det = detector)
}


MC_ocd_v6 <- function (Y, beta, sparsity, training_data = NA) {
  
  peak_stat <- future_map(1:length(Y), function(rep) {
    cat(rep, " ")
    ps <- c(0, 0, 0)
    y <- Y[[rep]]
    dim <- nrow(y)
    A <- matrix(0, dim, 1)
    
    if(is.na(training_data[1])){
      mu0 <- rep(0, dim)
    } else {
      mu0 <- apply(training_data[[rep]], 1, mean)
    }
    
    tail <- matrix(0, dim, floor(log2(dim)) * 2 + 4)
    for (i in 1:ncol(y)) {
      x_new <- y[, i] - mu0
      ret <- ocd_update(x_new, A, tail, beta, sparsity)
      A <- ret$A
      tail <- ret$tail
      ps <- pmax(ps, ret$stat)
    }
    return(ps)
  })
  
  peak_stat <- Reduce(rbind, peak_stat)
  colnames(peak_stat) <- c("diag", "off_d", "off_s")
  
  peak_stat
}
