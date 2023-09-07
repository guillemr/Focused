# devtools::install_github("abresler/nbastatR")
# Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2)
library(nbastatR)
library(tidyverse)
selectedSeasons <- c(2000:2022)
T_gamelog_reg <- suppressWarnings(game_logs(seasons = selectedSeasons, league = "NBA", result_types = "team", season_types = "Regular Season"))
View(head(T_gamelog_reg))

write_csv(T_gamelog_reg, file = "NBA_Regular_Season_2000_2022.csv")

