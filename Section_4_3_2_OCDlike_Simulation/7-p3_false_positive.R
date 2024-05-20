library(future)

CORES <- 10
plan(multicore, workers = CORES)
# Set up the cluster plan with the SSH session
# plan(cluster, workers = rep("SERVER ADDRESS", CORES))

# Script for the false positive rate with 3 dimentions

source("Section_4_3_2_OCDlike_Simulation/helper_functions.R")
source("MdFOCuS_R_implementation/MdFocus_MeanGaussian_md.R") ## code to test the main code

p <- 3
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
Y_train <- lapply(1:300, function(i) generate_sequence(n = 500, p = p, cp = 199, magnitude = 0, dens = 0, seed = 600 + i))

# data to train the monte-carlo treshold
Y_monte_carlo <- lapply(1:300, function(i) generate_sequence(n = N, p = p, cp = 500, magnitude = 0, dens = 0, seed = 123456 + i))


### pre-change mean known ###

##################################
###### FOCUS0 oracle #############
##################################

# tuning the threshold
focus0_mc <- future_map(Y_monte_carlo, function(y) {
  res <- FOCuS_multi_JMLR(y, c(Inf, Inf), mu0 = rep(0, p))
  data.frame(max = max(res$maxs), sum = max(res$sums))
}, .progress = T)
focus0_mc <- Reduce(rbind, focus0_mc)

alpha <- 0.01 / 2
thresholds$focus0 <- apply(focus0_mc, 2, quantile, probs = 1-alpha)


################################
##### md-focus0 oracle #########
################################

# tuning the threshold
md_focus0_mc <- future_map(Y_monte_carlo, function(y) {
  data <- t(y) # trasposing as the current prototype reads nxp (rather then pxn)
  res <- FocusCH(data, get_opt_cost = \(...) get_glo_opt(..., cost=lr_Focus0), threshold = Inf)
  - (res$opt.cost |> unlist())
}, .progress = T)
md_focus0_mc <- map_dbl(md_focus0_mc, max)

alpha <- 0.01
thresholds$mdfocus0 <- quantile(md_focus0_mc, probs = 1-alpha)

################################
####  md-focus0 part oracle ####
################################

# tuning the threshold
md_focus0_part_mc <- future_map(Y_monte_carlo, function(y) {
  data <- t(y) # trasposing as the current prototype reads nxp (rather then pxn)
  res <- FocusCH(data, get_opt_cost = \(...) get_3_opts_reg(..., cost=cost_lr_partial0), threshold = rep(Inf, 3))
  - (res$opt.cost |> reduce(rbind)) |> apply(2, max)
}, .progress = T)
md_focus0_part_mc <- reduce(md_focus0_part_mc, rbind)


alpha <- 0.01/3
thresholds$md_focus0_part <- apply(md_focus0_part_mc, 2, quantile, probs = 1-alpha)



############################
#####  ocd oracle ##########
############################

ocd_stat <- MC_ocd_v6(Y_monte_carlo, 1, "auto")

alpha <- 0.01/3
ocd_thres <- apply(ocd_stat, 2, quantile, prob = 1-alpha)

thresholds$ocd <- ocd_thres


#### pre-change mean unknown ####

##################################
###### FOCUS0 estimated ##########
##################################

# tuning the threshold
focus0est_mc <- future_map(1:REPS, function(i) {
  y <- Y_monte_carlo[[i]]
  y_tr <- Y_train[[i]]
  mu0hat <- rowMeans(y_tr)
  res <- FOCuS_multi_JMLR(y, c(Inf, Inf), mu0 = mu0hat)
  data.frame(max = max(res$maxs), sum = max(res$sums))
}, .progress = T)
focus0est_mc <- Reduce(rbind, focus0est_mc)

t_hat <- apply(focus0est_mc, 2, quantile, probs = .41)
t_multiplier <- cbind(focus0est_mc$max / t_hat["max"], focus0est_mc$sum / t_hat["sum"]) %>%
  apply(1, max) %>%
  quantile(probs = .41)
thresholds$focus0est <- t_hat * t_multiplier

# evaluating the empirical run length
focus0est_nc <- future_map(1:REPS, function(i) {
  y <- Y_nc[[i]]
  y_tr <- Y_train[[i]]
  mu0hat <- rowMeans(y_tr)
  res <- FOCuS_multi_JMLR(y, thresholds$focus0est, mu0 = mu0hat)
  res$t
}, .progress = T)
focus0est_nc <- focus0est_nc %>% unlist
focus0est_nc[focus0est_nc == -1] <- N
mean(focus0est_nc) # w\ current threshold 5028

####################
###### FOCUS  ######
####################

# tuning the threshold
focus_mc <- future_map(Y_monte_carlo, function(y) {
  res <- FOCuS_multi_JMLR(y, c(Inf, Inf))
  data.frame(max = max(res$maxs), sum = max(res$sums))
}, .progress = T)

focus_mc <- Reduce(rbind, focus_mc)
t_hat <- apply(focus_mc, 2, quantile, probs = .45)
t_multiplier <- cbind(focus_mc$max / t_hat["max"], focus_mc$sum / t_hat["sum"]) %>%
  apply(1, max) %>%
  quantile(probs = .45)
thresholds$focus <- t_hat * t_multiplier

# evaluating the empirical run length
focus_nc <- future_map(Y_nc, function(y) {
  res <- FOCuS_multi_JMLR(y, thresholds$focus)
  res$t
}, .progress = T)
focus_nc <- focus_nc %>% unlist
focus_nc[focus_nc == -1] <- N
mean(focus_nc) # w\ current 5086


########################
####  md-focus part ####
########################
# given that one penalty is included in the multiple penalty regimes,
# this simulation will just run the multiple one, and then 
# the result can be selected and tuned accordingly 

# tuning the threshold
md_focus_part_mc <- future_map(Y_monte_carlo, function(y) {
  data <- t(y) # trasposing as the current prototype reads nxp (rather then pxn)
  res <- FocusCH(data, get_opt_cost = get_partial_opt, threshold = rep(Inf, 3))
  - (res$opt.cost |> reduce(rbind)) |> apply(2, max)
}, .progress = T)
md_focus_part_mc <- reduce(md_focus_part_mc, rbind)


t_hat <- apply(md_focus_part_mc, 2, quantile, probs = .44)
t_multiplier <- as.data.frame(md_focus_part_mc) |> map2_df(t_hat, ~ .x / .y) %>%
  apply(1, max) %>%
  quantile(probs = .44)
thresholds$md_focus_part <- t_hat * t_multiplier

# evaluating the empirical run length
md_focus_part_nc <- future_map(Y_nc, function(y) {
  data <- t(y) # trasposing as the current prototype reads nxp (rather then pxn)
  res <- FocusCH(data, get_opt_cost = get_partial_opt, threshold = thresholds$md_focus_part)
  which(res$nb_at_step == 0)[1] - 1
}, .progress = T)
md_focus_part_nc <- md_focus_part_nc %>% map_dbl(~if_else(is.na(.x), N, .x[1]))
mean(md_focus_part_nc) # w\ current threshold 5024


########################
##### md-focus #########
########################

md_focus_mc <- md_focus_part_mc[, 3]
thresholds$mdfocus <- quantile(md_focus_mc, probs = .45)

# evaluating the empirical run length
md_focus_nc <- future_map(Y_nc, function(y) {
  data <- t(y) # trasposing as the current prototype reads nxp (rather then pxn)
  res <- FocusCH(data, get_opt_cost = \(...) get_glo_opt(..., cost=lr_Focus), threshold = thresholds$mdfocus)
  which(- (res$opt.cost |> unlist()) >= thresholds$mdfocus)
}, .progress = T)
md_focus_nc <- md_focus_nc %>% map_dbl(~if_else(is_empty(.x), N, .x[1]))
mean(md_focus_nc) # w\ current threshold 5037


###############################
#####  ocd estimated ##########
###############################

ocd_stat <- MC_ocd_v6(Y_monte_carlo, 1, "auto", training_data = Y_train)

avg_run_len <- 0

q_prob <- .68
while (avg_run_len < target_arl) {
  ocd_est_thres <- apply(ocd_stat, 2, quantile, prob = q_prob)
  
  res <- future_map(1:300, function(i) {
    y <- Y_nc[[i]]
    y_tr <- Y_train[[i]]
    
    ocd_det <- ocd_training(y_tr, ocd_est_thres)
    r <- ocd_detecting(y, ocd_det)
    r$t
  }, .progress = T)
  avg_run_len <- mean(unlist(res))
  print(avg_run_len)
  q_prob <- min(1, q_prob + 0.1)
}

thresholds$ocd_est <- ocd_est_thres # 6048.23



################################
###### saving the tresholds ####
################################

# remember to save the thresholds if you have updated anything!
save(thresholds, file = file)
