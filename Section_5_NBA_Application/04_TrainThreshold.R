#-----------------------------------------------------------------------------|
# SECTION 4.2.4 Simulating NBA Plus-Minus scores without changepoint          |
#                                                                             |
#-----------------------------------------------------------------------------|

################################################################################
##  Functions for Simulating the Plus-Minus score of matches                  ##
##  Pre-processing 2, Time Series Modeling                                    ##
################################################################################


# Description-------------------------------------------------------------|
# We generate :                                                           |
# (1) Skellam distribution without dependencies;                          |
# (2) Skellam distribution with dependencies (AR and sin);                |
# (3) the difference of two independent NegBin without dependencies;      |      
# (4) the difference of two independent NegBin with dependencies(AR and sin); 
# (5) Simulation of integer valued score uniformly at random between -80 and 80
# (6) Generation by re-sampling all past games before 2010                |
# (7) Generation by re-sampling past games of Cavaliers before 2010       |
#-------------------------------------------------------------------------|

# Remark------------------------------------------------------------------|
#                                                                         |
#-------------------------------------------------------------------------|
#parameters by default----------------------------------------------------|
set.seed(21) #@LUDA: I added this parameter
all_Mean <- list()
all_MeanAndVar <- list()

file_mean <- paste0("out_all_Mean_n_mc=", n_mc, ".RData")
file_var <- paste0("out_all_MeanAndVar_n_mc=", n_mc, ".RData")

save(all_Mean, file = file_mean)
save(all_MeanAndVar, file = file_var)

mu_skellams <- c(95, 105, 115)

# Generation with dependencies--------------------------------------------|
# ALL LIKE SKELLAM MEAN/VAR (GAUSS, SKELLAM AND MINUS NEG BIN)------------|
for (i in 1:length(sim_fn)) {   # add_fn and sim_fn generate in file "UTILS_FOR_TRAINING.R"
	for (j in 1:length(add_fn)) {
		for (k in 1:length(mu_skellams)) {
			suffix <- paste0(names(sim_fn)[i], "_", names(add_fn)[j], "_Mean=", mu_skellams[k])
			print(suffix) #@LUDA: remove this line?
			# generate time series with dependence
			sim.data <- function(n) sim_fn[[i]](n, mu_skellams[k], add_fn[[j]])

			# Attention : use parallel calculations mc.cores = mc.cores
			# For mean
			all_Mean[[suffix]] <- unlist(mclapply(rep(n, n_mc),
				FUN = sim_Max_ChangeInMean,
				sim.data = sim.data,
				mc.cores = mc.cores))
      # For mean/var
			all_MeanAndVar[[suffix]] <-  unlist(mclapply(rep(n, n_mc),
				FUN = sim_Max_ChangeInMeanAndVar,
  				sim.data = sim.data,
				min.size = min.size,
				minV = minV,
				mc.cores = mc.cores))
		}
	  #save in files
  	save(all_Mean, file = file_mean)
	  save(all_MeanAndVar, file = file_var)
	}
}

# Generation by uniform distribution-------------------------------------------|
#Simulation of integer valued score uniformly at random between -80 and 80
sim.data <- function(n) sample.int(160, size = n, replace = T)
suffix <- "Unif_Range_160"
# Attention : use parallel calculations mc.cores = mc.cores
# For mean
all_Mean[[suffix]] <- unlist(mclapply(rep(n, n_mc),
        FUN = sim_Max_ChangeInMean,
  			sim.data = sim.data,
  			mc.cores=6))
# For mean/var
all_MeanAndVar[[suffix]] <-  unlist(mclapply(rep(n, n_mc),
        FUN = sim_Max_ChangeInMeanAndVar,
  			sim.data = sim.data,
  			min.size = min.size,
  			minV = minV,
  			mc.cores = mc.cores))

# Generation by re-sampling past games-----------------------------------------|
### RESAMPLE FROM PAST NBA TEAM AND CAVALIERS

# read table with id of teams
dat_past <- dat %>% filter(yearSeason <= 2010)
past_all <- dat_past$plusminusTeam

dat_team_past <- dat_past %>% dplyr::filter(slugTeam == "CLE")
past_all_team  <- dat_team_past$plusminusTeam
#use all past games
sim.data <- function(n) sample(past_all, size = n, replace = T)
suffix <- "ResamplePast_All"
print(suffix) #remove?
# Attention : use parallel calculations mc.cores = mc.cores
# For mean
all_Mean[[suffix]] <- unlist(mclapply(rep(n, n_mc), 
                FUN = sim_Max_ChangeInMean,
  			        sim.data = sim.data, 
  			        mc.cores = mc.cores))
# For mean/var
all_MeanAndVar[[suffix]] <-  unlist(mclapply(rep(n, n_mc), 
        FUN = sim_Max_ChangeInMeanAndVar,
  			sim.data = sim.data, 
  			min.size = min.size, 
  			minV = minV, 
  			mc.cores = mc.cores))


#use all past games of Cavaliers
sim.data <- function(n) sample(past_all_team, size = n, replace = T)
suffix <- "ResamplePast_Cavaliers"
print(suffix) #remove?
# Attention : use parallel calculations mc.cores = mc.cores
# For mean
all_Mean[[suffix]] <- unlist(mclapply(rep(n, n_mc), 
        FUN = sim_Max_ChangeInMean,
  			sim.data = sim.data, 
  			mc.cores = mc.cores))
# For mean/var
all_MeanAndVar[[suffix]] <-  unlist(mclapply(rep(n, n_mc), 
        FUN = sim_Max_ChangeInMeanAndVar,
  			sim.data = sim.data, 
  			min.size = min.size, 
  			minV = minV, 
  			mc.cores = mc.cores))
#save

save(all_Mean, file = file_mean)
save(all_MeanAndVar, file = file_var)


################################################################################
########################### END ################################################
################################################################################

