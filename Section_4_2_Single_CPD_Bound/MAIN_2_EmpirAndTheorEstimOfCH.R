################################################################################
##     Figures Number of Facets and Vertices  of Convex Hull                  ##
##     Empirical results (using MAIN_1) and theoretical  estimations          ##
##     by Stirling numbers. Gaussian and Poisson Models                       ##
################################################################################

#packages-----------------------------------------------------------------|
library(ggplot2)
library(stats)
library(RColorBrewer)
library(ggpubr)
library(scales)

#function-----------------------------------------------------------------|
#function to calculate the partial sum of harmonic time series (for Stirling numbers)
part_sum_k_n <- function(k,n){
  data <- 1/(1: n[length(n)])^k
  res <-NULL
  for (i in 1:length(n))
    res <- c(res, sum(data[1:n[i]]))
  return (res)
}

#Gaussian Model###########################################################
#read results-------------------------------------------------------------|
results_fs <- list()
results_vs <- list()
#parameters---------------------------------------------------------------|
degree <- 10:23
n_ <-2^(degree)
dim_ <- 1 : 5
namedim_ <- paste("p =", dim_)
#number of simulations
nbSimus_ <- 100

u_n <- c(rep(2 * part_sum_k_n(1,n_), 
             each = nbSimus_),
         rep(2 * (part_sum_k_n(1,n_))^2 - 
               2 * part_sum_k_n(2,n_), 
             each = nbSimus_),
         rep(2 * (part_sum_k_n(1,n_))^3 - 
               6 * part_sum_k_n(1,n_) * part_sum_k_n(2,n_) +
               4 * part_sum_k_n(3,n_), 
             each = nbSimus_),
         rep(2*(part_sum_k_n(1,n_))^4 - 
               12 * (part_sum_k_n(1,n_))^2 * part_sum_k_n(2,n_) +
               16 * part_sum_k_n(1,n_) * part_sum_k_n(3,n_) +
               6 * (part_sum_k_n(2,n_))^2 -
               12 * part_sum_k_n(4,n_), 
             each = nbSimus_),
         rep(2 * (part_sum_k_n(1,n_))^5 -
               20 * (part_sum_k_n(1,n_))^3 * part_sum_k_n(2,n_) +
               40 * (part_sum_k_n(1,n_))^2 * part_sum_k_n(3,n_) -
               60 * part_sum_k_n(1,n_) * part_sum_k_n(4,n_) -
               40 * part_sum_k_n(2,n_) * part_sum_k_n(3,n_) +
               48 * part_sum_k_n(5,n_), 
             each = nbSimus_)
)

v_n <- c(rep(2 * part_sum_k_n(1,n_), each = nbSimus_),
         rep( (part_sum_k_n(1,n_))^2 - part_sum_k_n(2,n_) + 2, 
              each = nbSimus_),
         rep((1/3) * (part_sum_k_n(1,n_))^3 - part_sum_k_n(1,n_) * part_sum_k_n(2,n_) + 
               (2/3) * part_sum_k_n(3,n_) + 2 * part_sum_k_n(1,n_),
             each = nbSimus_),
         rep((1/12) * (part_sum_k_n(1,n_))^4 -
               (1/2) * (part_sum_k_n(1,n_))^2 * part_sum_k_n(2,n_) +
               (2/3) * part_sum_k_n(1,n_) * part_sum_k_n(3,n_) +
               (1/4) * (part_sum_k_n(2,n_))^2 -
               (1/2) * part_sum_k_n(4,n_) +
               (part_sum_k_n(1,n_))^2 -
               part_sum_k_n(2,n_) + 2, 
             each = nbSimus_),
         rep((1/60) * (part_sum_k_n(1,n_))^5 -
               (1/6) * (part_sum_k_n(1,n_))^3 * part_sum_k_n(2,n_) +
               (1/3) * (part_sum_k_n(1,n_))^2 * part_sum_k_n(3,n_) +
               (1/4) * part_sum_k_n(1,n_) * (part_sum_k_n(2,n_))^2 -
               (1/2) * part_sum_k_n(1,n_) * part_sum_k_n(4,n_) -
               (1/3) * part_sum_k_n(2,n_) * part_sum_k_n(3,n_) +
               (2/5) * part_sum_k_n(5,n_) +
               (1/3) * (part_sum_k_n(1,n_))^3 -
               part_sum_k_n(1,n_) * part_sum_k_n(2,n_) +
               (2/3) * part_sum_k_n(3,n_) + 
               2 * part_sum_k_n(1,n_), 
             each = nbSimus_)
)


for (j in  1 : length(dim_)) {
  results_fs[[j]] <- read.table(file = paste('Fs_dim',dim_[j],"gauss",nbSimus_,'.txt',sep = '_'), row.names = 1)
  results_vs[[j]] <- read.table(file = paste('Vs_dim',dim_[j],"gauss",nbSimus_,'.txt',sep = '_'), row.names = 1)
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


res_FS$u_n <-u_n
#plot
Plot_FS <- ggplot(res_FS, aes(x = n)) +
  geom_point(aes(y = y, color = Dimension), size = 1, shape = 1) +
  geom_line(aes(y = u_n, color = Dimension) , linetype = "dashed")  +
  scale_x_continuous("Number of data points, n-1", breaks=c(2^10,2^11,2^12,2^13, 2^14, 2^15, 2^16, 2^17,2^18,2^19,2^20,2^21,2^22,2^23),  trans = "log10",
                     labels = scales::math_format(2^.x, format = log2))+
  scale_y_continuous("Number of faces, Un", breaks=c(1,10,100,10^3, 10^4, 10^5, 10^6),  trans = "log10",
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
res_VS$v_n <-v_n

Plot_VS <- ggplot(res_VS, aes(x = n)) +
  geom_point(aes(y = y, color = Dimension), size = 1, shape = 1) +
  geom_line(aes(y = v_n, color = Dimension) , linetype = "dashed") +
  scale_x_continuous("Number of data points, n-1", breaks=c(2^10,2^11,2^12,2^13, 2^14, 2^15, 2^16, 2^17,2^18,2^19,2^20,2^21,2^22,2^23),  trans = "log10",
                     labels = scales::math_format(2^.x, format = log2))+
  scale_y_continuous("Number of vertices, Vn", breaks=c(1,10,100,10^3, 10^4),  trans = "log10",
                     labels = scales::math_format(10^.x, format = log10))+
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.position = "right") +
  scale_color_brewer(palette = 'Paired')

pdf(file = "Figure_Stirling_Gaussian_Model_Convex_hull_estimations.pdf",  width = 10, height = 4)
#par(mfrow=c(1,2))
#print(Plot_FS)
#print(Plot_VS)
ggarrange(Plot_FS, Plot_VS, ncol=2, common.legend = TRUE, legend = "right")
dev.off()


#Poisson Model#############################################################
#read results-------------------------------------------------------------|
results_fs <- list()
results_vs <- list()

for (j in  1 : length(dim_)) {
  results_fs[[j]] <- read.table(file = paste('Fs_dim',dim_[j],'poisson',nbSimus_,'.txt',sep = '_'), row.names = 1)
  results_vs[[j]] <- read.table(file = paste('Vs_dim',dim_[j],'poisson',nbSimus_,'.txt',sep = '_'), row.names = 1)
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
res_FS$u_n <-u_n
#plots
#plot
Plot_FS <- ggplot(res_FS, aes(x = n)) +
  geom_point(aes(y = y, color = Dimension), size = 1, shape = 1) +
  geom_line(aes(y = u_n, color = Dimension) , linetype = "dashed")  +
  scale_x_continuous("Number of data points, n-1", breaks=c(2^10,2^11,2^12,2^13, 2^14, 2^15, 2^16, 2^17,2^18,2^19,2^20,2^21,2^22,2^23),  trans = "log10",
                     labels = scales::math_format(2^.x, format = log2))+
  scale_y_continuous("Number of faces, Un", breaks=c(1,10,100,10^3, 10^4, 10^5, 10^6),  trans = "log10",
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


res_VS <- data.frame(n = data_frame_VS$n,
                     y = c(data_frame_VS[, 2],
                           data_frame_VS[, 3],
                           data_frame_VS[, 4],
                           data_frame_VS[, 5],
                           data_frame_VS[, 6]),
                     Dimension = rep(colnames(data_frame_VS[,-1]), each = nrow(data_frame_VS)))
res_VS$v_n <-v_n

Plot_VS <- ggplot(res_VS, aes(x = n)) +
  geom_point(aes(y = y, color = Dimension), size = 1, shape = 1) +
  geom_line(aes(y = v_n, color = Dimension) , linetype = "dashed") +
  scale_x_continuous("Number of data points, n-1", breaks=c(2^10,2^11,2^12,2^13, 2^14, 2^15, 2^16, 2^17,2^18,2^19,2^20,2^21,2^22,2^23),  trans = "log10",
                     labels = scales::math_format(2^.x, format = log2))+
  scale_y_continuous("Number of vertices, Vn", breaks=c(1,10,100,10^3, 10^4),  trans = "log10",
                     labels = scales::math_format(10^.x, format = log10))+
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.position = "right") +
  scale_color_brewer(palette = 'Paired')

pdf(file = "Figure_Stirling_Poisson_Model_Convex_hull_estimations.pdf",  width = 10, height = 4)
#par(mfrow=c(1,2))
#print(Plot_FS)
#print(Plot_VS)
ggarrange(Plot_FS, Plot_VS, ncol=2, common.legend = TRUE, legend = "right")
dev.off()
################################################################################
########################### END ################################################
################################################################################


