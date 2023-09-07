getDataTeam <- function(slugteam){
  dat_team <- dat %>% dplyr::filter(slugTeam == "CLE") %>%
  select(
    yearSeason,
    slugSeason,
    typeSeason,
    dateGame,
    nameTeam,
    slugTeam,
    isWin,
    ptsTeam,
    plusminusTeam
  ) %>%
  mutate(monthGame = format(as.Date(dateGame), "%Y-%m")) 
dat_team
}
