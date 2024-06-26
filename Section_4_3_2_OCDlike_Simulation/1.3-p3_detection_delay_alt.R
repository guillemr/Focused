library(future)
library(tidyr)
CORES <- 30
plan(multicore, workers = CORES)
# Set up the cluster plan with the SSH session
# plan(cluster, workers = rep("SERVER ADDRESS", CORES))

source("Section_4_3_2_OCDlike_Simulation/helper_functions.R")
source("MdFOCuS_R_implementation/MdFocus_MeanGaussian_md.R") ## code to test the main code


REPS <- 300

# setting parameters
p <- 3
N <- 5000
cp <- 1000
# loading thresholds from past simulations
load(paste0("Section_4_3_2_OCDlike_Simulation/results/thres", "p", p, "N", N, ".RData"))

# number of simulations with a change
sim_grid <- expand_grid(
  sim = 1:REPS,
  delta = c(2, 1, .75, .625, .5, 0.25, 0.125),  # magnitude of a change
  prop = 1:p/p,  # proportion of sequences with a change
  changepoint = cp,
  N = N
  )

file <- paste0("Section_4_3_2_OCDlike_Simulation/results/det_del_", "p", p, "N", N, "cp", cp, ".RData")
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

# FOCuS0 - oracle mean
focus0_res <- future_pmap(sim_grid, .f = function(delta, prop, changepoint, N, sim) {
  y <- generate_sequence(n = N, p = p, cp = changepoint, magnitude = delta, dens = prop, seed = sim)
  r <- FOCuS_multi_JMLR(y, thresholds$focus0, mu0 = rep(0, p))
  est <- ifelse(r$t == -1, N, r$t)
  data.frame(sim = sim, magnitude = delta, density = prop, algo = "FOCuS0", est = est, real = changepoint, N = N)
}, .progress = T)
focus0_res <- focus0_res %>% reduce(rbind)

runs_res$focus0 <- focus0_res

##################################
##### md-focus0 part oracle ######
##################################

# md-focus0 part oracle
md_focus0_part_res <- future_pmap(sim_grid, .f = function(delta, prop, changepoint, N, sim) {
  data <- t(generate_sequence(n = N, p = p, cp = changepoint, magnitude = delta, dens = prop, seed = sim)) # trasposing as the current prototype reads n x p (rather then p x n)
  res <- FocusCH(data, get_opt_cost = \(...) get_partial_opt(..., cost=cost_lr_partial0), threshold = thresholds$md_focus0_part)
  t <- which(res$nb_at_step == 0)[1]
  est <- if_else(is.na(t), N, t - 1)
  data.frame(sim = sim, magnitude = delta, density = prop, algo = "MdFOCuS0 MR", est = est, real = changepoint, N = N)
}, .progress = T)
md_focus0_part_res <- md_focus0_part_res %>% reduce(rbind)
runs_res$md_focus0_part <- md_focus0_part_res


############################
#####  ocd oracle ##########
############################

# ocd oracle
ocd_res <- future_pmap(sim_grid, .f = function(delta, prop, changepoint, N, sim) {
  y <- generate_sequence(n = N, p = p, cp = changepoint, magnitude = delta, dens = prop, seed = sim)
  ocd_det <- ocd_known(thresholds$ocd, rep(0, p), rep(1, p))
  r <- ocd_detecting(y, ocd_det)
  est <- r$t
  data.frame(sim = sim, magnitude = delta, density = prop, algo = "ocd (ora)", est = est, real = changepoint, N = N)
}, .progress = T)
ocd_res <- ocd_res %>% reduce(rbind)
runs_res$ocd <- ocd_res

###################### pre-change unkown ######################################

#####################
##### focus0 est ####
#####################

# 500 data points
focus0est_res <- future_pmap(sim_grid, .f = function(delta, prop, changepoint, N, sim) {
  y <- generate_sequence(n = N, p = p, cp = changepoint, magnitude = delta, dens = prop, seed = sim)
  y_tr <- generate_sequence(n = 500, p = p, cp = 199, magnitude = 0, dens = 0, seed = 600 + sim)
  mu0hat <- rowMeans(y_tr)
  r <- FOCuS_multi_JMLR(y, thresholds$focus0est, mu0 = mu0hat)
  est <- ifelse(r$t == -1, N, r$t)
  data.frame(sim = sim, magnitude = delta, density = prop, algo = "FOCuS0 (est 500)", est = est, real = changepoint, N = N)
}, .progress = T)
focus0est_res <- focus0est_res %>% reduce(rbind)
runs_res$focus0est <- focus0est_res

# 250 data points
focus0est_res <- future_pmap(sim_grid, .f = function(delta, prop, changepoint, N, sim) {
  y <- generate_sequence(n = N, p = p, cp = changepoint, magnitude = delta, dens = prop, seed = sim)
  y_tr <- generate_sequence(n = 250, p = p, cp = 199, magnitude = 0, dens = 0, seed = 600 + sim)
  mu0hat <- rowMeans(y_tr)
  r <- FOCuS_multi_JMLR(y, thresholds$focus0est_250, mu0 = mu0hat)
  est <- ifelse(r$t == -1, N, r$t)
  data.frame(sim = sim, magnitude = delta, density = prop, algo = "FOCuS0 (est 250)", est = est, real = changepoint, N = N)
}, .progress = T)
focus0est_res <- focus0est_res %>% reduce(rbind)
runs_res$focus0est_250 <- focus0est_res

# 100 data points
focus0est_res <- future_pmap(sim_grid, .f = function(delta, prop, changepoint, N, sim) {
  y <- generate_sequence(n = N, p = p, cp = changepoint, magnitude = delta, dens = prop, seed = sim)
  y_tr <- generate_sequence(n = 100, p = p, cp = 50, magnitude = 0, dens = 0, seed = 600 + sim)
  mu0hat <- rowMeans(y_tr)
  r <- FOCuS_multi_JMLR(y, thresholds$focus0est_100, mu0 = mu0hat)
  est <- ifelse(r$t == -1, N, r$t)
  data.frame(sim = sim, magnitude = delta, density = prop, algo = "FOCuS0 (est 100)", est = est, real = changepoint, N = N)
}, .progress = T)
focus0est_res <- focus0est_res %>% reduce(rbind)
runs_res$focus0est_100 <- focus0est_res


##################
##### focus ######
##################

# focus
focus_res <- future_pmap(sim_grid, .f = function(delta, prop, changepoint, N, sim) {
  y <- generate_sequence(n = N, p = p, cp = changepoint, magnitude = delta, dens = prop, seed = sim)
  r <- FOCuS_multi_JMLR(y, thresholds$focus)
  est <- ifelse(r$t == -1, N, r$t)
  data.frame(sim = sim, magnitude = delta, density = prop, algo = "FOCuS", est = est, real = changepoint, N = N)
}, .progress = T)
focus_res <- focus_res %>% reduce(rbind)
runs_res$focus <- focus_res


##################################
##### md-focus0 part oracle ######
##################################

md_focus_part_res <- future_pmap(sim_grid, .f = function(delta, prop, changepoint, N, sim) {
  data <- t(generate_sequence(n = N, p = p, cp = changepoint, magnitude = delta, dens = prop, seed = sim)) # trasposing as the current prototype reads n x p (rather then p x n)
  res <- FocusCH(data, get_opt_cost = get_partial_opt, threshold = thresholds$md_focus_part)
  t <- which(res$nb_at_step == 0)[1]
  est <- if_else(is.na(t), N, t - 1)
  data.frame(sim = sim, magnitude = delta, density = prop, algo = "MdFOCuS MR", est = est, real = changepoint, N = N)
}, .progress = T)
md_focus_part_res <- md_focus_part_res %>% reduce(rbind)
runs_res$md_focus_part <- md_focus_part_res



#######################
#####  ocd est ########
#######################


# ocd est
ocd_est_res <- future_pmap(sim_grid, .f = function(delta, prop, changepoint, N, sim) {
  y <- generate_sequence(n = N, p = p, cp = changepoint, magnitude = delta, dens = prop, seed = sim)
  y_tr <- generate_sequence(n = 500, p = p, cp = 199, magnitude = 0, dens = 0, seed = 600 + sim)
  ocd_det <- ocd_training(y_tr, thresholds$ocd_est)
  r <- ocd_detecting(y, ocd_det)
  est <- r$t
  data.frame(sim = sim, magnitude = delta, density = prop, algo = "ocd (est 500)", est = est, real = changepoint, N = N)
}, .progress = T)
ocd_est_res <- ocd_est_res %>% reduce(rbind)
runs_res$ocd_est <- ocd_est_res

# ocd est 250
ocd_est_res <- future_pmap(sim_grid, .f = function(delta, prop, changepoint, N, sim) {
  y <- generate_sequence(n = N, p = p, cp = changepoint, magnitude = delta, dens = prop, seed = sim)
  y_tr <- generate_sequence(n = 250, p = p, cp = 199, magnitude = 0, dens = 0, seed = 600 + sim)
  ocd_det <- ocd_training(y_tr, thresholds$ocd_est_250)
  r <- ocd_detecting(y, ocd_det)
  est <- r$t
  data.frame(sim = sim, magnitude = delta, density = prop, algo = "ocd (est 250)", est = est, real = changepoint, N = N)
}, .progress = T)
ocd_est_res <- ocd_est_res %>% reduce(rbind)
runs_res$ocd_est_250 <- ocd_est_res

# ocd est 100
ocd_est_res <- future_pmap(sim_grid, .f = function(delta, prop, changepoint, N, sim) {
  y <- generate_sequence(n = N, p = p, cp = changepoint, magnitude = delta, dens = prop, seed = sim)
  y_tr <- generate_sequence(n = 100, p = p, cp = 50, magnitude = 0, dens = 0, seed = 600 + sim)
  ocd_det <- ocd_training(y_tr, thresholds$ocd_est_100)
  r <- ocd_detecting(y, ocd_det)
  est <- r$t
  data.frame(sim = sim, magnitude = delta, density = prop, algo = "ocd (est 100)", est = est, real = changepoint, N = N)
}, .progress = T)
ocd_est_res <- ocd_est_res %>% reduce(rbind)
runs_res$ocd_est_100 <- ocd_est_res



### save past runs ###
save(runs_res, file = file)

### summaries ###
load(file)
tot_res <- Reduce(rbind, runs_res)

tot_summ <- tot_res %>% mutate(dd = est - cp, density = density * p) %>% 
  group_by(magnitude, density, algo) %>%
  summarise(avg_dd = mean(if_else(dd < 0, NA, dd), na.rm = T), fpr = mean(dd<0))
tot_summ %>% print(n=Inf)

wide_summary <- tot_summ %>% 
  pivot_wider(names_from = algo, values_from = avg_dd, id_cols = -fpr) 
View(wide_summary)

# unkown change
(wide_summary %>%
  select(magnitude, density, FOCuS0, "MdFOCuS0_part", "md-FOCuS0", "ocd (ora)")) 
# known change
 wide_summary %>%
  select(magnitude, density, "MdFOCuS MR", "ocd (est 250)", "ocd (est 100)")
