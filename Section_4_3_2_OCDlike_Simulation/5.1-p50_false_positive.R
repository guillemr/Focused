library(future)

CORES <- 30
plan(multicore, workers = CORES)
# Set up the cluster plan with the SSH session
# plan(cluster, workers = rep("romano@ma-res-romano.lancaster.ac.uk", 30))

# Script for the false positive rate with 100 dimensions

source("Section_4_3_2_OCDlike_Simulation/helper_functions.R")
source("MdFOCuS_R_implementation/MdFocus_MeanGaussian_md.R") ## code to test the main code

p <- 50
N <- 1e4
REPS <- 300
target_arl <- "fp_1"

file <- paste0("Section_4_3_2_OCDlike_Simulation/results/thres", "p", p, "N", target_arl, ".RData")
if (file.exists(file)) {
  # load the threshold list if the file is already there! 
  load(file)  
} else {
  # initialize null list for thresholds
  thresholds <- NULL
}

# data to estimate the mu0 value
Y_train <- future_map(1:300, function(i) generate_sequence(n = 500, p = p, cp = 199, magnitude = 0, dens = 0, seed = 600 + i))

# data to train the monte-carlo treshold
Y_monte_carlo <- future_map(1:300, function(i) generate_sequence(n = N, p = p, cp = 500, magnitude = 0, dens = 0, seed = 123456 + i))

### pre-change mean known ###

##################################
###### FOCUS0 oracle #############
##################################
if (T) {
  # tuning the threshold
  focus0_mc <- future_map(Y_monte_carlo, function(y) {
    res <- FOCuS_multi_JMLR(y, c(Inf, Inf), mu0 = rep(0, p))
    data.frame(max = max(res$maxs), sum = max(res$sums))
  }, .progress = T)
  focus0_mc <- Reduce(rbind, focus0_mc)
  
  alpha <- 0.01 / 2
  thresholds$focus0 <- apply(focus0_mc, 2, quantile, probs = 1-alpha)
}


################################
####  md-focus0 part oracle ####
################################

if (T) {
  # tuning the threshold
  md_focus0_part_mc <- future_map(Y_monte_carlo, function(y) {
    data <- t(y) # trasposing as the current prototype reads nxp (rather then pxn)
    res <- FocusCH_HighDim(data, get_opt_cost = \(...) get_partial_opt(..., cost=cost_lr_partial0, which_par = c(5, 25, 50)), threshold = rep(Inf, 5))
    - (res$opt.cost |> reduce(rbind)) |> apply(2, max)
  }, .progress = T)
  md_focus0_part_mc <- reduce(md_focus0_part_mc, rbind)
  
  save(md_focus0_part_mc, file="md_focus0_part1.RData")
  
} else {
  load("md_focus0_part1.RData")
}

if (T) {
  alpha <- 0.01/5
  thresholds$md_focus0_part <- apply(md_focus0_part_mc, 2, quantile, probs = 1-alpha)
}


###################################
####  md-focus0 1d part oracle ####
###################################

if (T) {
  # tuning the threshold
  md_focus0_1d_part_mc <- future_map(Y_monte_carlo, function(y) {
    data <- t(y) # trasposing as the current prototype reads nxp (rather then pxn)
    res <- FocusCH_HighDim(data, get_opt_cost = \(...) get_partial_opt(..., cost=cost_lr_partial0, which_par = c(5, 25, 50)), dim_indexes = as.list(1:ncol(data)), threshold = rep(Inf, 5))
    - (res$opt.cost |> reduce(rbind)) |> apply(2, max)
  }, .progress = T)
  md_focus0_1d_part_mc <- reduce(md_focus0_1d_part_mc, rbind)
  
  save(md_focus0_1d_part_mc, file="md_focus0_1d_part1.RData")
  
} else {
  load("md_focus0_1d_part1.RData")
}

if (T) {
  alpha <- 0.01/5
  thresholds$md_focus0_1d_part <- apply(md_focus0_1d_part_mc, 2, quantile, probs = 1-alpha)
}


# ############################
# #####  ocd oracle ##########
# ############################

# if (T) {
#   ocd_stat <- MC_ocd_v6(Y_monte_carlo, 1, "auto")
  
#   alpha <- 0.01/3
#   ocd_thres <- apply(ocd_stat, 2, quantile, prob = 1-alpha)

#   thresholds$ocd <- ocd_thres
  
# }


#### pre-change mean unknown ####

##################################
###### FOCUS0 estimated ##########
##################################
if (T) {
  # tuning the threshold
  focus0est_mc <- future_map(1:REPS, function(i) {
    y <- Y_monte_carlo[[i]]
    y_tr <- Y_train[[i]]
    mu0hat <- rowMeans(y_tr)
    res <- FOCuS_multi_JMLR(y, c(Inf, Inf), mu0 = mu0hat)
    data.frame(max = max(res$maxs), sum = max(res$sums))
  }, .progress = T)
  focus0est_mc <- Reduce(rbind, focus0est_mc)
  
  alpha <- 0.01 / 2
  thresholds$focus0est <- apply(focus0est_mc, 2, quantile, probs = 1-alpha)
  
}

####################
###### FOCUS  ######
####################
if (T) {
  # tuning the threshold
  focus_mc <- future_map(Y_monte_carlo, function(y) {
    res <- FOCuS_multi_JMLR(y, c(Inf, Inf))
    data.frame(max = max(res$maxs), sum = max(res$sums))
  }, .progress = T)
  focus_mc <- Reduce(rbind, focus_mc)

  alpha <- 0.01 / 2
  thresholds$focus <-  apply(focus_mc, 2, quantile, probs = 1-alpha)
  
}

########################
####  md-focus part ####
########################
# given that one penalty is included in the multiple penalty regimes,
# this simulation will just run the multiple one, and then 
# the result can be selected and tuned accordingly 

if(T) {# tuning the threshold
md_focus_part_mc <- future_map(Y_monte_carlo, function(y) {
  data <- t(y) # trasposing as the current prototype reads nxp (rather then pxn)
  res <- FocusCH_HighDim(data, get_opt_cost = \(...) get_partial_opt(..., which_par = c(5, 25, 50)), threshold = rep(Inf, 5))
  - (res$opt.cost |> reduce(rbind)) |> apply(2, max)
}, .progress = T)
md_focus_part_mc <- reduce(md_focus_part_mc, rbind)

save(md_focus_part_mc, file="md_focus_part1.RData")
} else {
  load("md_focus_part1.RData")
}

if(T){

  alpha <- 0.01/5
  thresholds$md_focus_part <- apply(md_focus_part_mc, 2, quantile, probs = 1-alpha)

}


###################################
####  md-focus 1d part oracle ####
###################################

if (T) {
  # tuning the threshold
  md_focus_1d_part_mc <- future_map(Y_monte_carlo, function(y) {
    data <- t(y) # trasposing as the current prototype reads nxp (rather then pxn)
    res <- FocusCH_HighDim(data, get_opt_cost = \(...) get_partial_opt(..., which_par = c(5, 25, 50)), dim_indexes = as.list(1:ncol(data)), threshold = rep(Inf, 5))
    - (res$opt.cost |> reduce(rbind)) |> apply(2, max)
  }, .progress = T)
  md_focus_1d_part_mc <- reduce(md_focus_1d_part_mc, rbind)
  
  save(md_focus_1d_part_mc, file="md_focus_1d_part1.RData")
  
} else {
  load("md_focus_1d_part1.RData")
}


if (T) {
  alpha <- 0.01/5
  thresholds$md_focus_1d_part <- apply(md_focus_1d_part_mc, 2, quantile, probs = 1-alpha)
}


# ###############################
# #####  ocd estimated ##########
# ###############################
# if(T) {
# ocd_stat <- MC_ocd_v6(Y_monte_carlo, 1, "auto", training_data = Y_train)


# alpha <- 0.01 / 3
# thresholds$ocd_est <-  apply(ocd_stat, 2, quantile, probs = 1 - alpha)

################################
###### saving the tresholds ####
################################

# remember to save the thresholds if you have updated anything!
save(thresholds, file = file)
