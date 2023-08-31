#-----------------------------------------------------------------------------|
# SECTION 4.2.4 Simulating NBA Plus-Minus scores without changepoint          |
#                                                                             |
#-----------------------------------------------------------------------------|

################################################################################
##  NBA Application. Pre-processing                                           ##
################################################################################


#functions---------------------------------------------------------------|
#------------------------------------------------------------------------|
#' @title transform_diff
#'
#' @description difference between score of the team and score of the opponent 
#' @param score_team score of the team
#' @param score_oppo score of the opponent
#' @return difference between score of the team and score of the opponent 
transform_diff <- function(score_team, score_oppo) {
  score_team - score_oppo
}

#------------------------------------------------------------------------|
#' @title transform_logdiff
#'
#' @description the difference between the logarithm of the team's score and the logarithm of the opponent's score
#' @param score_team score of the team
#' @param score_oppo score of the opponent
#' @return the difference between the logarithm of the team's score and the logarithm of the opponent's score
transform_logdiff <- function(score_team, score_oppo) {
  log(score_team) - log(score_oppo)
}

#------------------------------------------------------------------------|
#' @title getData
#'
#' @description Function to return information of all matches of team in time range
#' @param nickname nickname of the team
#' @param dat date of game
#' @param first_season first season of this team
#' @param last_season last season of this team
#' @param transform  type of score transformation (by default, transform_diff)
#' @return Function to return information of all matches of team in time range
getData <- function(nickname, dat, first_season, last_season, transform = transform_diff) {
  #pre-processing
  id_team <- teams$TEAM_ID[teams$NICKNAME == nickname] 
  home_matches <- which(id_team == dat$TEAM_ID_home)
  away_matches <- which(id_team == dat$TEAM_ID_away)

  all_matches_home <- data.frame(SEASON = dat$SEASON[home_matches], 
                                 YEAR = dat$YEAR[home_matches], 
                                 MONTH = dat$MONTH[home_matches], 
                                 DAY = dat$DAY[home_matches], 
                                 score_oppo = dat$PTS_away[home_matches], 
                                 score_team = dat$PTS_home[home_matches], 
                                 home = TRUE)

  all_matches_away <- data.frame(SEASON = dat$SEASON[away_matches], 
                                 YEAR = dat$YEAR[away_matches], 
                                 MONTH = dat$MONTH[away_matches], 
                                 DAY = dat$DAY[away_matches], 
                                 score_oppo = dat$PTS_home[away_matches], 
                                 score_team = dat$PTS_away[away_matches], 
                                 home = FALSE)

  #bind and order
  dat_matches <- rbind(all_matches_away, all_matches_home)
  dat_matches <- dat_matches[order(dat_matches$YEAR, dat_matches$MONTH, dat_matches$DAY), ]

  dat_matches <- dat_matches[dat_matches$SEASON >= first_season & dat_matches$SEASON <= last_season, ]
  dat_matches$s_all  <- transform(dat_matches$score_team, dat_matches$score_oppo)
  dat_matches$s_home <- dat_matches$s_all
  dat_matches$s_away <- dat_matches$s_all
  dat_matches$s_home[!dat_matches$home] <- NA
  dat_matches$s_away[ dat_matches$home] <- NA
  
  return(dat_matches)
}

#------------------------------------------------------------------------|
#' @title winner 
#'
#' @description the difference between the logarithm of the team's score and the logarithm of the opponent's score
#' @param score score
#' @param m normalization parameter
#' @return return TRUE if score > 0
winner <- function(score, m = NA) {
  score > 0
}

#------------------------------------------------------------------------|
#' @title winner_step
#'
#' @description 
#' @param score score
#' @param m normalization parameter
#' @return 
winner_step <- function(score, m = NA) {
   out <- rep(1, length(score))
   out[is.na(score)] <- 0
   return(out)
}

#------------------------------------------------------------------------|
#' @title difference
#'
#' @description 
#' @param score score
#' @param m normalization parameter (by default, m = 0.494)
#' @return normalization of score (for second e-detectors)
difference <- function(score, m = 0.494) {
  (score + 80)/ (160 * m) - 1
}

#------------------------------------------------------------------------|
#' @title difference_step
#'
#' @description 
#' @param score score
#' @param m normalization parameter (by default, m = 0.494)
#' @return 
difference_step <- function(score, m = 0.494) {
  out <- ((score + 80)/ (160 * m) - 1)^2 
  out[is.na(score)] <- 0
  return(out)
}

#------------------------------------------------------------------------|
#' @title getCand_1d 
#'
#' @description Function to return the candidate set by the Quickhull algorithm 
#' for Gaussian Model (change in mean)
#' @param dat_matches 
#' @param score_function 
#' @param step_function 
#' @param m normalization parameter (by default, m = 0.494)
#' @return the candidate set by the Quickhull algorithm for Gaussian Model (change in mean)
getCand_1d <- function(dat_matches, score_function, step_function, m = 0.494) {
  #prepare for qhull (get random walk)
  dat_analysis <- data.frame(match = step_function(dat_matches$s_all, m), 
                             score = score_function(dat_matches$s_all, m))
  walk_analysis <- apply(dat_analysis, 2, cumsum)[-nrow(dat_analysis), ]
  #Convex hull analysis
  hull_analysis <- convhulln(walk_analysis[, 1:2])   # get convex hull
  cpt_cand <- sort(unique(as.vector(hull_analysis))) # get vertices
  return(cpt_cand)
}

#------------------------------------------------------------------------|
#' @title getCand_meanvar
#'
#' @description Function to return the candidate set by the Quickhull algorithm 
#' for Gaussian Model
#' @param dat_matches 
#' @param onlymean if TRUE - change in mean, else change in mean and variance
#' @return the candidate set by the Quickhull algorithm for Gaussian Model 
getCand_meanvar <- function(dat_matches, onlymean = T) {
  #prepare for qhull (get random walk)
  dat_analysis <- data.frame(match = 1, score = (dat_matches$s_all), score2 = (dat_matches$s_all)^2)
  walk_analysis <- apply(dat_analysis, 2, cumsum)[-nrow(dat_analysis), ]
  #Convex hull analysis
  if(onlymean)  hull_analysis <- convhulln(walk_analysis[, 1:2]) # get convex hull (change in mean)
  else hull_analysis <- convhulln(walk_analysis[, 1:3]) # get convex hull (change in mean and var)
  cpt_cand <- sort(unique(as.vector(hull_analysis))) # get vertices
  return(cpt_cand)
}

#------------------------------------------------------------------------|
#' @title getCand_2d 
#'
#' @description Function to return the candidate set by the Quickhull algorithm 
#' @param dat_matches 
#' @param score_function 
#' @param step_function 
#' @param m normalization parameter (by default, m = 0.494)
#' @return the candidate set by the Quickhull algorithm for Gaussian Model (change in mean)
getCand_2d <- function(dat_matches, score_function, step_function, m = 0.494) {
  #prepare for qhull (get random walk)
  dat_analysis <- data.frame(home = step_function(dat_matches$s_home, m), 
                             away = step_function(dat_matches$s_away, m), 
                             scoreHome = score_function(dat_matches$s_home, m), 
                             scoreAway = score_function(dat_matches$s_away, m))
  dat_analysis$scoreHome[!dat_matches$home] <- 0.
  dat_analysis$scoreAway[ dat_matches$home] <- 0.
  walk_analysis <- apply(dat_analysis, 2, cumsum)[-nrow(dat_analysis), ]
  #Convex hull analysis
  hull_analysis <- convhulln(walk_analysis[, 1:4]) # get convex hull 
  cpt_cand <- sort(unique(as.vector(hull_analysis))) # get vertices
  return(cpt_cand)
}
#------------------------------------------------------------------------|
#' @title getEst
#'
#' @description  
#' @param dat_matches 
#' @param tau 
#' @return 
getEst <- function(dat_matches, tau) {
  dat_matches$segment <- 2
  dat_matches$segment[1:tau] <- 1
  colnames(ave_win)[-1] <- gsub("s_", "w_", colnames(ave_win)[-1])
  ave_score <- aggregate(dat_matches[, 8:10], by = list(dat_matches$segment), FUN = mean, na.rm = T)
  ave_win <- aggregate(dat_matches[, 8:10] > 0, by = list(dat_matches$segment), FUN = mean, na.rm = T)
  colnames(ave_win)[-1] <- gsub("s_", "w_", colnames(ave_win)[-1])
  cbind(signif(ave_score[, -1], 2), signif(ave_win[, -1], 2))
}

#------------------------------------------------------------------------|
#' @title getEst
#'
#' @description  
#' @param dat_matches 
#' @param tau 
#' @return 
getData_spec <- function(dat_matches, sigma = 1) {
  dat_analysis <- data.frame(home = dat_matches$home + 0, 
                             away = (!dat_matches$home) + 0, 
                             scoreHome = dat_matches$s_home/sigma, 
                             scoreAway = dat_matches$s_away/sigma)
  dat_analysis$scoreHome[!dat_matches$home] <- 0.
  dat_analysis$scoreAway[ dat_matches$home] <- 0.
  return(dat_analysis)
}



################################################################################
########################### END ################################################
################################################################################


