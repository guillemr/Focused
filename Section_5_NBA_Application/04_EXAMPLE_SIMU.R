################################################################################
##                    Figure 8 (Right) Modeling Dependencies. Examples             ##
################################################################################

#packages-----------------------------------------------------------------|
library(extraDistr)
library(ggpubr) 
library(RColorBrewer)
library(ggplot2)
#source-----------------------------------------------------------------|
source("04_Utils_For_Training.R")
#parameters-----------------------------------------------------------------|
#parameters-----------------------------------------------------------------|
n <- 1000
mu_skellams <- c(95, 105, 115)

#loop
for (i in 1 : length(sim_fn)) {   # add_fn and sim_fn generate in file "UTILS_FOR_TRAINING.R"
  for (j in 1 : length(add_fn)) {
    for (k in 1 : length(mu_skellams)) {
      add_dep <- add_fn[[j]]
      add_to_mu1 <- add_dep(n)
      mu_skellam <- mu_skellams[k]
      x <- rskellam(n, mu1 = mu_skellam + add_to_mu1, mu2 = mu_skellam)
      
      data <- data.frame(match_nb = 1:n ,
                         PTS_Differential = x,
                         Game = c("Lost", "Won")[(x > 0) + 1],
                         mean = add_to_mu1)
      
      name <- paste0(names(sim_fn)[i], "-", names(add_fn)[j], "-mu=", mu_skellam)
      p <- ggplot(data, aes(x = match_nb, y = PTS_Differential, fill = Game)) +
        geom_bar(stat = "identity", position = "identity") +
        scale_fill_manual(values = c("#D04A61","#1F78B4"))+ 
        geom_line(aes(x = match_nb, y = mean), color = "black", linetype = "dashed")+
       # labs(title = paste0(names(sim_fn)[i], "-", names(add_fn)[j], "-mu=", mu_skellam))+
        scale_y_continuous(limits = c(-60, 60)) + 
        ylab("Point differential") +
        xlab("Number of matches") +
        theme_bw() +
        theme(text = element_text(size = 13),
              legend.title = element_text(size = 13),
              legend.text = element_text(size = 13),
              legend.position = "right") 
      pdf(file = paste0("figure/simu_example/", name, ".pdf"),  width = 5.5, height = 4) 
      print(p)
      dev.off()
      #ggsave(paste0("figure/simu_example/", name, ".pdf"))
    }}}

