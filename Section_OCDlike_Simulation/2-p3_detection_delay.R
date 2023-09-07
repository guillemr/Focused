source("simulations/helper_functions.R")
source("md_focus0_draft/simpleCH.R") ## code with the main function
source("md_focus0_draft/simpleCH_utils.R") ## code to test the main code

CORES <- 16
REPS <- 300

# setting parameters
p <- 3
N <- 5000
# loading thresholds from past simulations
load(paste0("simulations/thres", "p", p, "N", N, ".RData"))

# number of simulations with a change
sim_grid <- expand.grid(
  delta = c(2, 1, .5, 0.25),  # magnitude of a change
  prop = 1:p/p,  # proportion of sequences with a change
  changepoint = 500,
  N = N
)

file <- paste0("simulations/det_del_", "p", p, "N", N, ".RData")
if (file.exists(file)) {
  # load the existing simulation
  load(file)  
} else {
  ### or initialize a list for saving the individuals runs ###
  runs_res <- NULL
}


##########################
##### focus0 oracle ######
##########################

focus0_res <- NULL

for (s in 1:nrow(sim_grid)) {
  simu <- sim_grid[s, ]
  Y <- mclapply(1:REPS, function(i) generate_sequence(n = simu$N, p = p, cp = simu$changepoint, magnitude = simu$delta, dens = simu$prop, seed = i),
                mc.cores = CORES)
  
  # FOCuS0 - oracle mean
  res <- mclapply(1:REPS, function(i) {
    y <- Y[[i]]
    r <- FOCuS_multi_JMLR(y, thresholds$focus0, mu0 = rep(0, p))
    ifelse(r$t == -1, simu$N, r$t)
  }, mc.cores = CORES)
  res <- unlist(res)
  focus0_res[[s]] <- data.frame(sim = 1:REPS, magnitude = simu$delta, density = simu$prop, algo = "FOCuS0", est = res, real = simu$changepoint, N = simu$N)
}
focus0_res <- focus0_res %>% reduce(rbind)
runs_res$focus0 <- focus0_res


#############################
##### md-focus0 oracle ######
#############################

md_focus0_res <- NULL

for (s in 1:nrow(sim_grid)) {
  simu <- sim_grid[s, ]
  print(simu)
  Y <- mclapply(1:REPS, function(i) generate_sequence(n = simu$N, p = p, cp = simu$changepoint, magnitude = simu$delta, dens = simu$prop, seed = i),
                mc.cores = CORES)
  
  # FOCuS0 - oracle mean
  res <- mclapply(1:REPS, function(i) {
    data <- t(Y[[i]]) # trasposing as the current prototype reads n x p (rather then p x n)
    res <- FocusCH(data, fun.cost=lr_Focus0, common_difference_step=1, common_ratio_step=2, first_step_qhull=ncol(data)+2, threshold = thresholds$mdfocus0)
    over_the_t <- which(-res$opt.cost >= thresholds$mdfocus0)
    if_else(is_empty(over_the_t), simu$N, over_the_t[1])
  }, mc.cores = CORES)
  res <- unlist(res)
  md_focus0_res[[s]] <- data.frame(sim = 1:REPS, magnitude = simu$delta, density = simu$prop, algo = "md-FOCuS0", est = res, real = simu$changepoint, N = simu$N)
}
md_focus0_res <- md_focus0_res %>% reduce(rbind)
runs_res$md_focus0 <- md_focus0_res


############################
#####  ocd oracle ##########
############################

ocd_res <- NULL
for (s in 1:nrow(sim_grid)) {
  simu <- sim_grid[s, ]
  print(simu)
  Y <- mclapply(1:REPS, function(i) generate_sequence(n = simu$N, p = p, cp = simu$changepoint, magnitude = simu$delta, dens = simu$prop, seed = i),
                mc.cores = CORES)
  
  # ocd - oracle mean
  res <- mclapply(1:REPS, function(i) {
    y <- Y[[i]]
    ocd_det <- ocd_known(thresholds$ocd, rep(0, p), rep(1, p))
    r <- ocd_detecting(y, ocd_det)
    r$t
  }, mc.cores = CORES)
  res <- unlist(res)
  ocd_res[[s]] <- data.frame(sim = 1:REPS, magnitude = simu$delta, density = simu$prop, algo = "ocd (ora)", est = res, real = simu$changepoint, N = simu$N)
}
runs_res$ocd <- ocd_res %>% reduce(rbind)

###################### pre-change unkown ######################################

#####################
##### focus0 est ####
#####################

focus0est_res <- NULL

for (s in 1:nrow(sim_grid)) {
  simu <- sim_grid[s, ]
  Y <- mclapply(1:REPS, function(i) generate_sequence(n = simu$N, p = p, cp = simu$changepoint, magnitude = simu$delta, dens = simu$prop, seed = i),
                mc.cores = CORES)
  
  # FOCuS0 - estimated mean
  res <- mclapply(1:REPS, function(i) {
    y <- Y[[i]]
    y_tr <- Y_train[[i]]
    mu0hat <- rowMeans(y_tr)
    r <- FOCuS_multi_JMLR(y, thresholds$focus0est, mu0 = mu0hat)
    ifelse(r$t == -1, simu$N, r$t)
  }, mc.cores = CORES)
  res <- unlist(res)
  focus0est_res[[s]] <- data.frame(sim = 1:REPS, magnitude = simu$delta, density = simu$prop, algo = "FOCuS0 (est)", est = res, real = simu$changepoint, N = simu$N)
}
focus0est_res <- focus0est_res %>% reduce(rbind)
runs_res$focus0est <- focus0est_res



##################
##### focus ######
##################

focus_res <- NULL

for (s in 1:nrow(sim_grid)) {
  simu <- sim_grid[s, ]
  Y <- mclapply(1:REPS, function(i) generate_sequence(n = simu$N, p = p, cp = simu$changepoint, magnitude = simu$delta, dens = simu$prop, seed = i),
                mc.cores = CORES)
  
  # focus - oracle mean
  res <- mclapply(1:REPS, function(i) {
    y <- Y[[i]]
    r <- FOCuS_multi_JMLR(y, thresholds$focus)
    ifelse(r$t == -1, simu$N, r$t)
  }, mc.cores = CORES)
  res <- unlist(res)
  focus_res[[s]] <- data.frame(sim = 1:REPS, magnitude = simu$delta, density = simu$prop, algo = "FOCuS", est = res, real = simu$changepoint, N = simu$N)
}
focus_res <- focus_res %>% reduce(rbind)
runs_res$focus <- focus_res


####################
##### md-focus######
####################

md_focus_res <- NULL

for (s in 1:nrow(sim_grid)) {
  simu <- sim_grid[s, ]
  print(simu)
  Y <- mclapply(1:REPS, function(i) generate_sequence(n = simu$N, p = p, cp = simu$changepoint, magnitude = simu$delta, dens = simu$prop, seed = i),
                mc.cores = CORES)
  
  # focus - oracle mean
  res <- mclapply(1:REPS, function(i) {
    data <- t(Y[[i]]) # trasposing as the current prototype reads n x p (rather then p x n)
    res <- FocusCH(data, fun.cost=lr_Focus, common_difference_step=1, common_ratio_step=2, first_step_qhull=ncol(data)+2, threshold = thresholds$mdfocus)
    over_the_t <- which(-res$opt.cost >= thresholds$mdfocus)
    if_else(is_empty(over_the_t), simu$N, over_the_t[1])
  }, mc.cores = CORES)
  res <- unlist(res)
  md_focus_res[[s]] <- data.frame(sim = 1:REPS, magnitude = simu$delta, density = simu$prop, algo = "md-FOCuS", est = res, real = simu$changepoint, N = simu$N)
}
md_focus_res <- md_focus_res %>% reduce(rbind)
runs_res$md_focus <- md_focus_res


#######################
#####  ocd est ########
#######################

ocd_est_res <- NULL
for (s in 1:nrow(sim_grid)) {
  simu <- sim_grid[s, ]
  print(simu)
  Y <- mclapply(1:REPS, function(i) generate_sequence(n = simu$N, p = p, cp = simu$changepoint, magnitude = simu$delta, dens = simu$prop, seed = i),
                mc.cores = CORES)
  Y_train <- lapply(1:300, function(i) generate_sequence(n = 500, p = p, cp = 199, magnitude = 0, dens = 0, seed = 600 + i))
  
  
  # ocd - est mean
  res <- mclapply(1:REPS, function(i) {
    y <- Y[[i]]
    y_tr <- Y_train[[i]]
    ocd_det <- ocd_training(y_tr, thresholds$ocd_est)
    r <- ocd_detecting(y, ocd_det)
    r$t
  }, mc.cores = CORES)
  res <- unlist(res)
  ocd_est_res[[s]] <- data.frame(sim = 1:REPS, magnitude = simu$delta, density = simu$prop, algo = "ocd (est)", est = res, real = simu$changepoint, N = simu$N)
}
runs_res$ocd_est <- ocd_est_res %>% reduce(rbind)



### save past runs ###
save(runs_res, file = file)

### summaries ###
load(file)
tot_res <- Reduce(rbind, runs_res)

tot_summ <- tot_res %>% mutate(dd = est - 500, density = density * p) %>% 
  group_by(magnitude, density, algo) %>%
  summarise(avg_dd = mean(if_else(dd < 0, NA, dd), na.rm = T), fpr = mean(dd<0))
tot_summ %>% print(n=Inf)

tot_summ %>% pivot_wider(names_from = algo, values_from = avg_dd, id_cols = -fpr)
