#-----------------------------------------------------------------------------|
# SECTION 4.2 Illustration on NBA Plus-Minus                                  |
#         4.2.4 Simulating NBA Plus-Minus scores without changepoint          |
#-----------------------------------------------------------------------------|

################################################################################
##                    Figure 6 Evolution of the CUMSUM Score                  ##
## Gaussian Model(change in mean) and Gaussian Model(change in mean and var)  ##                
################################################################################
rm(list = ls())
#packages-----------------------------------------------------------------|
library(geometry)
library(parallel)
library(ggplot2)
library(extraDistr)
library(base)
library(RColorBrewer)
library(scales)
library(ggpubr) # for plot
library(tidyverse)
#include files------------------------------------------------------------|
# ATTENTION: check the path to the directory
source("../Md_Focus_R_Functions/MdFocus_MeanVarGaussian_1d.R")               # MdFOCuS for Gaussian Model(change in mean and variance)
source("../Md_Focus_R_Functions/MdFocus_MeanGaussian_md.R")                # MdFOCuS for Gaussian Model(change in mean)
source("04_Utils_For_Training.R")      # Pre-processing 1: Modeling dependencies and functions for time series generation
source("04_Read_NBA.R")

#parameters by default----------------------------------------------------|
n         <- 10^3     # number of data points
n_mc      <- 10^3
min.size  <- 2
minV      <- 1 
## variance of home matches per team overall season much larger than 1




#Get all matches and time---------------------------------------------------|
dat <- read.table("NBA_Regular_Season_2000_2022.csv", sep = ",", header = T)

#Variance of plus minus ---------------------------------------------------|
#Get an estimate of the variance for all teams and all seasons
all_variances_per_season <- unlist(lapply(unique(dat$slugTeam), FUN=function(slugTeam){
  dat_team <- getDataTeam(slugTeam)
  tapply(dat_team$plusminusTeam, INDEX=dat_team$yearSeason, FUN=var)
}
))


summary(all_variances_per_season)
##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  97.02  146.36  168.06  170.68  187.05  247.38 

#Mean of plus minus ---------------------------------------------------|
#Get an estimate of the mean plusminus per season

  
summary( tapply(dat$ptsTeam, INDEX=(dat$slugSeason), FUN=mean))
##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  93.40   97.10   99.92  101.07  104.13  112.09 

# read table of games for cavaliers
dat_team <- getDataTeam("CLE")
dat_team <- dat_team %>% filter(yearSeason > 2010 & yearSeason <= 2018)


## consecutive equal plus minus score  
sum(diff(dat_team$plusminusTeam)==0)
which(diff(dat_team$plusminusTeam)==0)
sum(dat_team$plusminusTeam[1:638] == dat_team$plusminusTeam[2:639] & dat_team$plusminusTeam[1:638] == dat_team$plusminusTeam[3:640])

## season change
season_change <- which(diff(dat_team$yearSeason)!=0) # number of matches

#Transformation of score for random walk --------------------------------| 
x <- dat_team$plusminusTeam
data1  <- cbind(x, x^2)          

#Detect change by MdFOCuS for Gaussian model (change in mean)------------| 
res_Mean <- FocusCH(data1[, 1, drop = F], 
                    fun.cost = lr_Focus,
                    common_difference_step = 1, 
                    common_ratio_step = 2, 
                    first_step_qhull = ncol(data1) + 5)

#Detect change by  by MdFOCuS for Gaussian model (change in mean and variance)-------| 
res_MeanAndVar <- FocusCH_var(data1, 
                              fun.cost = lr_Focus_var, 
                              min.size = min.size, 
                              minV = minV, 
                              common_difference_step = 1, 
                              common_ratio_step = 2, 
                              first_step_qhull = ncol(data1) + 5)


# Data frame with the colons: number of matches, optimal cost (cumsum) by Gaussian model(change in mean) 
# and optimal cost (cumsum) by Gaussian model(change in mean and variance)
data2 <- data.frame(nb_match = 1:nrow(dat_team),
                    cumsum_mean = -res_Mean$opt.cost, 
                    cumsum_meanandvar = -res_MeanAndVar$opt.cost)
colnames(data2) <- c('Number of matches', 'Change in mean', 'Change in mean and variance')

# for mean "#6A3D9A"; for mean and  "#E31A1C"

#get plot for  Gaussian model(change in mean)
p_mean <- ggplot(data2, aes(x = data2[, 1], y = data2[, 2])) + 
  geom_line(color = "#6A3D9A") +  
  ylab("Cumsum Statistics") +
  xlab("Number of matches") +
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.position = "right") +
  geom_vline(xintercept = 335.5, linetype = "dotted")

#save plot
pdf(file = "figure/Figure_cumsum_change_in_mean.pdf",  width = 5, height = 4) 
print(p_mean)
dev.off()


#get plot for  Gaussian model(change in mean and variance)
p_meanvar <- ggplot(data2, aes(x = data2[, 1], y = data2[, 3])) + 
  geom_line(color =  "#E31A1C") +  
  ylab("Cumsum Statistics") +
  xlab("Number of matches") +
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.position = "right") +
  geom_vline(xintercept = 335.5, linetype = "dotted")
#save plot
pdf(file = "figure/Figure_cumsum_change_in_mean_var.pdf",  width = 5, height = 4) 
print(p_meanvar)
dev.off() 

#save two plots
p_mean_meanvar <- ggarrange(p_mean, p_meanvar, ncol = 2)
pdf(file = "figure/Figure_Cumsum.pdf",  width = 10, height = 4) 
print(p_mean_meanvar)
dev.off() 



################################################################################
##TEST: Simulating NBA Plus-Minus scores without changepoint                  ##
## Gaussian Model(change in mean) and Gaussian Model(change in mean and var)  ## 
## Threshold and detection time                                               ##
################################################################################

# Test description--------------------------------------------------------|
# Idea: To explore a few natural optimistic and pessimistic options to  
# model the Plus-Minus score  and to check how this affect the threshold
# of the two models ('change in mean' and 'change in mean and variance') 
# and in the time at which they detect the change.
#-------------------------------------------------------------------------|

# Remark------------------------------------------------------------------|
# We consider the following simulatons:                                   
# (1) Skellam distribution without dependencies;                          
# (2) Skellam distribution with dependencies (AR and sin);                
# (3) the difference of two independent NegBin without dependencies;            
# (4) the difference of two independent NegBin with dependencies(AR and sin); 
# (5) Simulation of integer valued score uniformly at random between -80 and 80
# (6) Generation by re-sampling all past games before 2010                
# (7) Generation by re-sampling past games of Cavaliers before 2010       
# See details in Pre-processing 1 (UTILS_FOR_TRAINING.R) and 
# Pre-processing 2 (Train_MeanVar.R)
#-------------------------------------------------------------------------|

# Get training data for analysis------------------------------------------|
#(Gaussian model (change in mean) and Gaussian model (change in mean and variance) 
file_mean <- paste0("out_all_Mean_n_mc=", n_mc, ".RData")
file_var  <- paste0("out_all_MeanAndVar_n_mc=", n_mc, ".RData")
mc.cores  <-6
if(!file.exists(file_mean )) source("04_TrainThreshold.R")    # Generate training (attention : 3h with mc.cores = 6)   
load(file = file_mean)
load(file = file_var)

# Recover 95% thresholds-------------------------------------------------|
thrs_Mean <- sapply(all_Mean, function(x) quantile(x, prob = 1 - 0.05))               # for mean model 
thrs_MeanAndVar <- sapply(all_MeanAndVar, function(x) quantile(x, prob = 1 - 0.05))   # for mean/var model 

data <- data.frame(
  value = c(thrs_Mean, thrs_MeanAndVar),
  group = rep(c("Change in mean", "Change in mean and variance"), each = length(thrs_Mean))
)
################################################################################
##                    Figure 7 (Left) Threshold obtain for                    ##
## Gaussian Model(change in mean) and   Gaussian Model(change in mean and var)##                
################################################################################
# Creating the histogram
p_thrs <- ggplot(data, aes(x = value, fill = group)) +
  geom_histogram(binwidth = 0.1, position = "identity", alpha = 0.7) +
  scale_x_continuous("Frequency", breaks=c(10^0,10^1,10^2,10^3,10^4, 10^5),  trans = "log10",
                     labels = scales::math_format(10^.x, format = log10))+
  scale_fill_discrete(name = "Type of model") +
  ylab("Value of treshold") +
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 13),
        legend.position = "right") +
  scale_color_brewer(palette = 'Paired') +
  labs(title = "Histogram (in log-scale)") + 
  scale_fill_manual(values = c("#6A3D9A","#E31A1C")) +
  labs(fill = "Type of Gaussian Model") 

#save plot
pdf(file = "figure/Figure_Treshold_by_training.pdf",  width = 10, height = 4) 
print(p_thrs)
dev.off()


################################################################################
##                    Figure 7 (Right) Detection time obtain for              ##
## Gaussian Model(change in mean) and   Gaussian Model(change in mean and var)##                
################################################################################
time_Mean <- sapply(thrs_Mean, FUN = function(thrs) min(c(which(thrs < -(res_Mean$opt.cost)), 10^3)))
time_MeanAndVar <- sapply(thrs_MeanAndVar, FUN = function(thrs) min(c(which(thrs < -(res_MeanAndVar$opt.cost)), 10^3)))

range(time_Mean)
## [1]  330 1000
mean(time_Mean==1000)## = not detection
range(time_MeanAndVar)
## [1] 339 368

### percentage of simulation models for which the detecting is later with the mean model
mean(time_Mean > time_MeanAndVar)
### 0.942029

data3 <- data.frame(
  value = c(time_Mean, time_MeanAndVar),
  group = rep(c("Change in mean", "Change in mean and variance"), each = length(thrs_Mean))
)
#get plot
p_time_detect <- ggplot(data3, aes(x = value, fill = group)) +
  geom_histogram(binwidth = 5, position = "identity", alpha = 0.7) +
  scale_fill_discrete(name = "Type of Gaussian Model") +
  labs(x = "Value of time detection", y = "Frequency",title = "Histogram") +
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 13),
        legend.position = "right") + 
  scale_fill_manual(values = c("#6A3D9A","#E31A1C")) +
  labs(fill = "Type of Gaussian Model") +
    geom_vline(xintercept = 854, linetype = "dotted")
#save plot
pdf(file = "figure/Figure_Time_detection_NBA_Application.pdf",  width = 10, height = 4) 
print(p_time_detect)
dev.off()

p_thrs_time_detect <- ggarrange(p_thrs , p_time_detect, ncol = 2,common.legend = TRUE, legend = "bottom")
#save two plots
pdf(file = "figure/Figure_Treshold_and_Time_detection_by_training.pdf",  width = 10, height = 4) 
print(p_thrs_time_detect)
dev.off() 

##############################################################
## Detected change                                          ##
##############################################################
change_Mean <- sapply(time_Mean, FUN = function(t) res_Mean$opt.change[t])
change_MeanAndVar<- sapply(time_MeanAndVar, FUN = function(t) res_MeanAndVar$opt.change[t])

## always detect a change at 279 a bit after the start of the season at 312
range(change_Mean, na.rm = T) ## 279
unique(res_Mean$opt.change[320:nrow(dat_team)]) ## 279 

range(change_MeanAndVar, na.rm = T) ## 279
unique(res_MeanAndVar$opt.change[320:nrow(dat_team)]) ## 279 

season_change <- which(diff(dat_team$yearSeason)!=0)
season_change[4] ## change between "2013-14" "2014-15" match 312

#------------------------------------------------------------------------------|

# Transformation of  the data as in e-detectors first to get convex-hull points
#Convex-hull for e-detectors
### Win-rate
x <- dat_team$plusminusTeam
data_ <- cbind(1, x > 0 - 0.49)
cs_data_ <- apply(data_, 2, cumsum)
nb_point <- sapply(10:nrow(cs_data_), FUN = function(n) length(unique(as.vector(convhulln(cs_data_[1:n, ])))))
range(nb_point) ## 5 to 16 (compare to the number of lambda values used in e-detectors)


sort(unique(as.vector(convhulln(cs_data_))))
# [1]   1   4   7  53  70 176 246 279 300 324 351 626 639 640

# Transformation of  the data as in e-detectors second to get convex-hull points
m <- 0.494
x <- dat_team$plusminusTeam
x <- (x + 80)/160
x <- x/m - 1
data_ <- cbind(1, x, x^2) ## the one is for the mean-var 2:3 only for
cs_data_ <- apply(data_, 2, cumsum)
nb_point <- sapply(10:nrow(cs_data_), FUN = function(n) length(unique(as.vector(convhulln(cs_data_[1:n, 2:3])))))
range(nb_point) ## 5 to 18 (compare to the number of lambda values used in e-detectors)
sort(unique(as.vector(convhulln(cs_data_[, 2:3]))))
## closest to 312 : 300 and 351
## [1]   1   4  37  65  70 176 246 279 300 351 624 626 632 639 64

## Mean and Variance model
nb_point <- sapply(10:nrow(cs_data_), FUN = function(n) length(unique(as.vector(convhulln(cs_data_[1:n, ])))))
range(nb_point) ## 5 to 55
sort(unique(as.vector(convhulln(cs_data_))))
## closest to 312 : 300 and 319
## [1]   1   2   4   7   8   9  10  11  15  16  37  40  48  51  53  65  67  69  70
## [20]  72 150 158 176 184 195 246 264 266 279 280 300 304 319 350 351 353 356 392
## [39] 396 554 566 595 597 624 625 626 632 633 635 639 640


closest_to_312 <- sapply(10:nrow(cs_data_), FUN = function(n) min(abs(312 -unique(as.vector(convhulln(cs_data_[1:n, ])))) ))
which(closest_to_312 == 0) + 9

## all points on the hull for the mean and e-detector model should be on the
## hull of the mean and variance model
meanAndVar_Cand <- sort(unique(as.vector(convhulln(cs_data_))))
mean_Cand <- sort(unique(as.vector(convhulln(cs_data_[, 1:2]))))
eDetect_Cand <- sort(unique(as.vector(convhulln(cs_data_[, 2:3]))))
eDetect_Cand %in% meanAndVar_Cand
mean_Cand %in% meanAndVar_Cand



### Rather than the arrival of Lebron James does match a change of General Manager
### David Griffin Wikie : On February 6, 2014, he was named the acting general manager for the Cavaliers, after the firing of Chris Grant
### see https://en.wikipedia.org/wiki/David_Griffin_(basketball)
### Fairly good end of season from win better point difference and win-rate
dat_team[279:280, ]
#   yearSeason slugSeason     typeSeason   dateGame            nameTeam
# 279       2014    2013-14 Regular Season 2014-02-05 Cleveland Cavaliers
# 280       2014    2013-14 Regular Season 2014-02-07 Cleveland Cavaliers
x <- dat_team$plusminusTeam
before_change <- x[230:279]
between_change_and_newseason <- x[280:312]

mean(before_change)
mean(between_change_and_newseason)
t.test(before_change, between_change_and_newseason)

################################################################################
##  Figure 6: Plus-Minus score of the Cavaliers from season 2013 to 2014 as a bar plot. 
## Losses are in light blue and wins are in dark blue. 
##Ends of seasons are represented with vertical dashed lines(grey) and seasons are shown on the bottom with grey arrows.              ##
################################################################################
x_ <- dat_team$plusminusTeam
match_nb <- (season_change[3] + 1): season_change[5]
data <- data.frame(
  match_nb = match_nb ,
  PTS_Differential = x_[match_nb],
  Game = c("Won", "Lost")[(x_[match_nb] > 0) + 1]
)
#get plot
arrow_change <- 279+.5
p <- ggplot(data, aes(x = match_nb, y = PTS_Differential, fill = Game)) +
  geom_bar(stat = "identity", position = "identity") +
  scale_fill_manual(values = c("#A6CEE3", "#1F78B4")) +
  labs(
    #title = "Point differential",
    x = "Matches ordered chronogicaly",
    y = "Point differential"
) +
  geom_vline(xintercept = c(arrow_change ), color = "#E31A1C", linetype = "dashed") +
  geom_vline(xintercept = c(season_change[4]+0.5), color = "grey", linetype = "dashed") +
  annotate(
    geom = "segment",
    x = arrow_change,
    xend = arrow_change,
    y = 27,  # Extend the arrow from bottom to top of the plot area
    yend = 10,
    color = "#E31A1C",
    arrow = arrow(type = "closed", length = unit(0.10, "inches"))
  ) +
  annotate(
    geom = "segment",
    x = season_change[3]+0.5,
    xend = season_change[4],
    y = -50,  # Extend the arrow from bottom to top of the plot area
    yend = -50,
    color = "grey",
    linewidth = 1.5,
    arrow = arrow(type = "closed", length = unit(0.1, "inches"))
  ) +
  annotate(
    geom = "segment",
    x = season_change[4]+0.5,
    xend = season_change[5],
    y = -50,  # Extend the arrow from bottom to top of the plot area
    yend = -50,
    color = "grey",
    linewidth = 1.5,
    arrow = arrow(type = "closed", length = unit(0.1, "inches"))
  ) +
  theme_bw()+
  geom_text(aes(x = arrow_change, y = 35, label = "Change detected  \n between   \n the 5th and 7th of Febuary"), size = 4) +
  geom_text(aes(x = mean(season_change[3:4]), y = -47, label = "Season 2013-14"), size = 4, color ="grey") +
  geom_text(aes(x = mean(season_change[4:5]), y = -47, label = "Season 2014-15"), size = 4, color = "grey") + 
  theme(text = element_text(size = 13),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 13),
        legend.position = "right")
#save plot
pdf(file = "figure/Figure_point_differential_2013_2014.pdf",  width = 10, height = 4) 
print(p)
dev.off()

t.test(x_[(season_change[3]+1):279], x_[280:season_change[4]])
## p-value of 0.016


################################################################################
## Figure 4: Plus-Minus score of the Cavaliers from season 2010 to season 2018 as a bar plot. 
## Losses are in light blue and wins are in dark blue.
## Ends of seasons are represented with vertical dashed lines and seasons are shown on the bottom with grey arrows.
################################################################################
x_ <- dat_team$plusminusTeam
match_nb <- 1:length(x_)
data <- data.frame(
  match_nb = match_nb ,
  PTS_Differential = x_[match_nb],
  Game = c("Won", "Lost")[(x_[match_nb] > 0)+1]
)

season_change <- which(diff(dat_team$yearSeason)!=0) # number of matches
season_change <- c(0, season_change, length(x_))
p_all <- ggplot(data, aes(x = match_nb, y = PTS_Differential, fill = Game)) +
  geom_bar(stat = "identity", position = "identity") +
  scale_fill_manual(values = c("#A6CEE3", "#1F78B4")) +
  labs(
    #title = "Point differential",
    x = "Matches ordered chronogicaly",
    y = "Point differential"
) + theme_bw()

for(i in  1:(length(season_change)-1)){
p_all <- p_all + annotate(
    geom = "segment",
    x = season_change[i]+0.5,
    xend = season_change[i+1],
    y = -50,  # Extend the arrow from bottom to top of the plot area
    yend = -50,
    color = c("grey", "darkgrey")[i %% 2 + 1],
    linewidth=1,
    arrow = arrow(type = "closed", length = unit(0.1, "inches"))
  ) +
   geom_vline(xintercept = c(season_change[i]+0.5), color = "grey", linetype = "dashed")
}


p_all <- p_all +  geom_text(aes(x = mean(season_change[1:2]), y = -47, label = paste0("", 2009+1, "-", 9+2)) , color ="grey")+
geom_text(aes(x = mean(season_change[1:2+1]), y = -47, label = paste0("", 2009+2, "-", 9+3)) , color ="grey")+
geom_text(aes(x = mean(season_change[1:2+2]), y = -47, label = paste0("", 2009+3, "-", 9+4)) , color ="grey")+
geom_text(aes(x = mean(season_change[1:2+3]), y = -47, label = paste0("", 2009+4, "-", 9+5)) , color ="grey")+
 geom_text(aes(x = mean(season_change[1:2+4]), y = -47, label = paste0("", 2009+5, "-", 9+6)) , color ="grey")+
 geom_text(aes(x = mean(season_change[1:2+5]), y = -47, label = paste0("", 2009+6, "-", 9+7)) , color ="grey")+
 geom_text(aes(x = mean(season_change[1:2+6]), y = -47, label = paste0("", 2009+7, "-", 9+8)) , color ="grey")+
 geom_text(aes(x = mean(season_change[1:2+7]), y = -47, label = paste0("", 2009+8, "-", 9+9)) , color ="grey")+
 geom_text(aes(x = mean(season_change[1:2+8]), y = -47, label = paste0("", 2009+9, "-", 9+10)) , color ="grey") +
  theme(text = element_text(size = 13),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 13),
        legend.position = "right")
#save plot
pdf(file = "figure/Figure_Seasons_2010_2018.pdf",  width = 10, height = 4) 
print(p_all)
dev.off()




################################################################################
########################### END ################################################
################################################################################




