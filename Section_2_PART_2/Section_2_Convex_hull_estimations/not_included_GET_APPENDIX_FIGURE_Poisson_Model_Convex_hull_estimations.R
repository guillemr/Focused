################################################################################
##     Figure (not included) Number of Facets and Vertices  of Convex Hull    ##
##     Empirical results and theoretical asymptotic estimations               ##
##                               Poisson model                                ##
################################################################################

#packages-----------------------------------------------------------------|
library(ggplot2)
library(stats)
library(RColorBrewer)
library(ggpubr)
library(scales)

#parameters (by default)--------------------------------------------------|
dim_ <- 1:5
nbSimus_ <- 100
namedim_ <- c ("p = 1","p = 2","p = 3","p = 4","p = 5")

#read results-------------------------------------------------------------|
results_fs <- list()
results_vs <- list()
for (j in  1 : length(dim_)) {
  results_fs[[j]] <- read.table(file = paste('Fs_dim',
                                             dim_[j],
                                             'poisson',
                                             nbSimus_,
                                             '.txt',
                                             sep = '_'), row.names = 1)

  results_vs[[j]] <- read.table(file = paste('Vs_dim',
                                             dim_[j],
                                             'poisson',
                                             nbSimus_,
                                             '.txt',
                                             sep = '_'), row.names = 1)
}

data_frame_FS <- data.frame(n = results_fs[[1]][1], matrix(NA, dim(results_fs[[1]][1])[1], length(dim_)))
colnames(data_frame_FS) <- c('n', namedim_)

data_frame_VS <- data.frame(n = results_vs[[1]][1], matrix(NA, dim(results_vs[[1]][1])[1], length(dim_)))
colnames(data_frame_VS) <- c('n', namedim_)

for (i in 1:length(dim_) ) {
  data_frame_FS[1:dim(results_fs[[i]])[1],i+1] <- results_fs[[i]][2]
  data_frame_VS[1:dim(results_fs[[i]])[1],i+1] <-results_vs[[i]][2]
}

res_FS <- data.frame(n = data_frame_FS$n,
                     y = c(data_frame_FS[, 2],
                           data_frame_FS[, 3],
                           data_frame_FS[, 4],
                           data_frame_FS[, 5],
                           data_frame_FS[, 6]),
                     Dimension = rep(colnames(data_frame_FS[,-1]), each = nrow(data_frame_FS)))
#brewer.pal(n=5,"Set1")
#[1] "#E41A1C" "#377EB8" "#4DAF4A" "#984EA3" "#FF7F00"

data_frame_FS <- data.frame(n = results_fs[[1]][1], matrix(NA, dim(results_fs[[1]][1])[1], length(dim_)))
colnames(data_frame_FS) <- c('n', namedim_)

data_frame_VS <- data.frame(n = results_vs[[1]][1], matrix(NA, dim(results_vs[[1]][1])[1], length(dim_)))
colnames(data_frame_VS) <- c('n', namedim_)

for (i in 1:length(dim_) ) {
  data_frame_FS[1:dim(results_fs[[i]])[1],i+1] <- results_fs[[i]][2]
  data_frame_VS[1:dim(results_fs[[i]])[1],i+1] <-results_vs[[i]][2]
}

res_FS <- data.frame(n = data_frame_FS$n,
                     y = c(data_frame_FS[, 2],
                           data_frame_FS[, 3],
                           data_frame_FS[, 4],
                           data_frame_FS[, 5],
                           data_frame_FS[, 6]),
                     Dimension = rep(colnames(data_frame_FS[,-1]), each = nrow(data_frame_FS)))
#brewer.pal(n=5,'Paired')
#"#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99"

Plot_FS <- ggplot(res_FS, aes(x = n, y = y, color = Dimension)) +
  geom_point( size = 1, shape = 1) +
  geom_line(aes(y = 2*log(n)), linetype = "dashed",color = "#A6CEE3") +
  geom_line(aes(y = 2*(log(n))^2), linetype = "dashed", color = "#1F78B4") +
  geom_line(aes(y = 2*(log(n))^3), linetype = "dashed", color = "#B2DF8A") +
  geom_line(aes(y = 2*(log(n))^4), linetype = "dashed", color = "#33A02C") +
  geom_line(aes(y = 2*(log(n))^5), linetype = "dashed", color = "#FB9A99") +
  scale_x_continuous("Number of data points, n-1", breaks=c(2^10,2^11,2^12,2^13, 2^14, 2^15, 2^16, 2^17,2^18,2^19,2^20,2^21,2^22,2^23),  trans = "log10",
                     labels = scales::math_format(2^.x, format = log2))+
  scale_y_continuous("Number of facets, u\u2099", breaks=c(1,10,100,10^3, 10^4, 10^5, 10^6),  trans = "log10",
                     labels = scales::math_format(10^.x, format = log10)) +
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.position = "right") +
  scale_color_brewer(palette = 'Paired')

res_VS <- data.frame(n = data_frame_VS$n,
                     y = c(data_frame_VS[, 2],
                           data_frame_VS[, 3],
                           data_frame_VS[, 4],
                           data_frame_VS[, 5],
                           data_frame_VS[, 6]),
                     Dimension = rep(colnames(data_frame_VS[,-1]), each = nrow(data_frame_VS)))

Plot_VS <- ggplot(res_VS, aes(x = n, y = y, color = Dimension)) +
  geom_point(size = 1, shape = 1) +
  geom_line(aes(y = (2/factorial(1))*log(n)),  linetype = "dashed",color= "#A6CEE3") +
  geom_line(aes(y = (2/factorial(2))*(log(n))^2),  linetype = "dashed", color = "#1F78B4") +
  geom_line(aes(y = (2/factorial(3))*(log(n))^3),  linetype = "dashed", color = "#B2DF8A") +
  geom_line(aes(y = (2/factorial(4))*(log(n))^4), linetype = "dashed", color = "#33A02C") +
  geom_line(aes(y = (2/factorial(5))*(log(n))^5),  linetype = "dashed", color = "#FB9A99") +
  scale_x_continuous("Number of data points, n-1", breaks=c(2^10,2^11,2^12,2^13, 2^14, 2^15, 2^16, 2^17,2^18,2^19,2^20,2^21,2^22,2^23),  trans = "log10",
                     labels = scales::math_format(2^.x, format = log2))+
  scale_y_continuous("Number of vertices, v\u2099", breaks=c(1,10,100,10^3, 10^4),  trans = "log10",
                     labels = scales::math_format(10^.x, format = log10))+
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.position = "right") +
  scale_color_brewer(palette = 'Paired')
#print(Plot_FS)
#print(Plot_VS)

pdf(file = "Figure_Convex_hull_asymptotic_estimations_Poisson(not included in article).pdf",  width = 8.5, height = 4)
ggarrange(Plot_FS, Plot_VS, ncol=2, common.legend = TRUE, legend = "right")
dev.off()

