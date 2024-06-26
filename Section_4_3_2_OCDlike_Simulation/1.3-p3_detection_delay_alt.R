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
  select(magnitude, density, "MdFOCuS MR", "ocd (est 500)", "ocd (est 100)")
