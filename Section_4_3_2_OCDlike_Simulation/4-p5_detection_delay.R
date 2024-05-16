source("Section_4_3_2_OCDlike_Simulation/helper_functions.R")
source("MdFOCuS_R_implementation/MdFocus_MeanGaussian_md.R") ## code to test the main code

#source("md_focus0_draft/simpleCH.R") ## code with the main function
#source("md_focus0_draft/simpleCH_utils.R") ## code to test the main code

CORES <- 10
plan(multicore, workers = CORES)

REPS <- 300

# setting parameters
p <- 5
N <- 5000
# loading thresholds from past simulations
load(paste0("Section_4_3_2_OCDlike_Simulation/results/thres", "p", p, "N", N, ".RData"))

# number of simulations with a change
sim_grid <- expand_grid(
  sim = 1:REPS,
  delta = c(2, 1, .5, 0.25),  # magnitude of a change
  prop = 1:p/p,  # proportion of sequences with a change
  changepoint = 500,
  N = N
  )

file <- paste0("Section_4_3_2_OCDlike_Simulation/results/det_del_", "p", p, "N", N, ".RData")
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
  data.frame(sim = sim, magnitude = delta, density = prop, algo = "FOCuS0_alt", est = est, real = changepoint, N = N)
}, .progress = T)
focus0_res <- focus0_res %>% reduce(rbind)

runs_res$focus0 <- focus0_res



#############################
##### md-focus0 oracle ######
#############################

# md-focus0 oracle
md_focus0_res <- future_pmap(sim_grid, .f = function(delta, prop, changepoint, N, sim) {
  data <- t(generate_sequence(n = N, p = p, cp = changepoint, magnitude = delta, dens = prop, seed = sim)) # trasposing as the current prototype reads n x p (rather then p x n)
  res <- FocusCH(data, get_opt_cost = \(...) get_glo_opt(..., cost=lr_Focus0), threshold = thresholds$mdfocus0)
  t <- which(res$nb_at_step == 0)[1]
  est <- if_else(is.na(t), N, t - 1)
  data.frame(sim = sim, magnitude = delta, density = prop, algo = "md-FOCuS0", est = est, real = changepoint, N = N)
}, .progress = T)
md_focus0_res <- md_focus0_res %>% reduce(rbind)
runs_res$md_focus0 <- md_focus0_res


##################################
##### md-focus0 part oracle ######
##################################

# md-focus0 part oracle
md_focus0_part_res <- future_pmap(sim_grid, .f = function(delta, prop, changepoint, N, sim) {
  data <- t(generate_sequence(n = N, p = p, cp = changepoint, magnitude = delta, dens = prop, seed = sim)) # trasposing as the current prototype reads n x p (rather then p x n)
  res <- FocusCH(data, get_opt_cost = \(...) get_partial_opt(..., cost=cost_lr_partial0), threshold = thresholds$md_focus0_part)
  t <- which(res$nb_at_step == 0)[1]
  est <- if_else(is.na(t), N, t - 1)
  data.frame(sim = sim, magnitude = delta, density = prop, algo = "MdFOCuS0_part", est = est, real = changepoint, N = N)
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

focus0est_res <- future_pmap(sim_grid, .f = function(delta, prop, changepoint, N, sim) {
  y <- generate_sequence(n = N, p = p, cp = changepoint, magnitude = delta, dens = prop, seed = sim)
  y_tr <- generate_sequence(n = 500, p = p, cp = 199, magnitude = 0, dens = 0, seed = 600 + sim)
  mu0hat <- rowMeans(y_tr)
  r <- FOCuS_multi_JMLR(y, thresholds$focus0est, mu0 = mu0hat)
  est <- ifelse(r$t == -1, N, r$t)
  data.frame(sim = sim, magnitude = delta, density = prop, algo = "FOCuS0 (est)", est = est, real = changepoint, N = N)
}, .progress = T)
focus0est_res <- focus0est_res %>% reduce(rbind)
runs_res$focus0est <- focus0est_res


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


####################
##### md-focus######
####################

# md-focus
md_focus_res <- future_pmap(sim_grid, .f = function(delta, prop, changepoint, N, sim) {
  data <- t(generate_sequence(n = N, p = p, cp = changepoint, magnitude = delta, dens = prop, seed = sim)) # trasposing as the current prototype reads n x p (rather then p x n)
  res <- FocusCH(data, get_opt_cost = \(...) get_glo_opt(..., cost=lr_Focus), threshold = thresholds$mdfocus)
  t <- which(res$nb_at_step == 0)[1]
  est <- if_else(is.na(t), N, t - 1)
  data.frame(sim = sim, magnitude = delta, density = prop, algo = "md-FOCuS", est = est, real = changepoint, N = N)
}, .progress = T)
md_focus_res <- md_focus_res %>% reduce(rbind)
runs_res$md_focus <- md_focus_res



##################################
##### md-focus0 part oracle ######
##################################

md_focus_part_res <- future_pmap(sim_grid, .f = function(delta, prop, changepoint, N, sim) {
  data <- t(generate_sequence(n = N, p = p, cp = changepoint, magnitude = delta, dens = prop, seed = sim)) # trasposing as the current prototype reads n x p (rather then p x n)
  res <- FocusCH(data, get_opt_cost = get_partial_opt, threshold = thresholds$md_focus_part)
  t <- which(res$nb_at_step == 0)[1]
  est <- if_else(is.na(t), N, t - 1)
  data.frame(sim = sim, magnitude = delta, density = prop, algo = "MdFOCuS_part", est = est, real = changepoint, N = N)
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
  data.frame(sim = sim, magnitude = delta, density = prop, algo = "ocd (est)", est = est, real = changepoint, N = N)
}, .progress = T)
ocd_est_res <- ocd_est_res %>% reduce(rbind)
runs_res$ocd_est <- ocd_est_res


### save past runs ###
save(runs_res, file = file)

### summaries ###
load(file)
tot_res <- Reduce(rbind, runs_res)

tot_summ <- tot_res %>% mutate(dd = est - 500, density = density * p) %>% 
  group_by(magnitude, density, algo) %>%
  summarise(avg_dd = mean(if_else(dd < 0, NA, dd), na.rm = T), fpr = mean(dd<0))
tot_summ %>% print(n=Inf)

wide_summary <- tot_summ %>% 
  pivot_wider(names_from = algo, values_from = avg_dd, id_cols = -fpr) 

# unkown change
(wide_summary %>%
  select(magnitude, density, FOCuS0, "MdFOCuS0_part", "md-FOCuS0", "ocd (ora)")) 
# known change
wide_summary %>%
  select(magnitude, density, FOCuS, "FOCuS0 (est)", "MdFOCuS_part", "md-FOCuS", "ocd (est)")
