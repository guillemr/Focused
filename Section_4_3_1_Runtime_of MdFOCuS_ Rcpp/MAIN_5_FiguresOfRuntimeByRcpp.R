################################################################################
##     Figure Empirical Run time as the function of alpha                     ##                             
##                      Gaussian Model (MdFOCuS0, p=2 and 3)                  ##
################################################################################
# Remark------------------------------------------------------------------|
# We use the results obtained in MAIN_3 and MAIN_4                        |
#-------------------------------------------------------------------------|

#packages-----------------------------------------------------------------|
library(ggplot2)
library(stats)
library(RColorBrewer)
library(ggpubr)
library(scales)
#parameters---------------------------------------------------------------|
alpha_ <- c(1+2^(-4:0), 3, 2^(2:7))
N <-10^5
NbSimus <- 100
Method <-c("FOCuS0")
Cost <- c("gauss")
#data frame---------------------------------------------------------------|
table_alpha <- data.frame(alpha = c(alpha_, alpha_, alpha_),
                          Runtime = rep(NA, 3*length(alpha_)), Dimension = c(rep("p = 1", length(alpha_)), rep("p = 2", length(alpha_)), rep("p = 3", length(alpha_))))
colnames(table_alpha) <- c('\u03b1', 'Run time, seconds', 'Dimension')
table_alpha [1:length(alpha_), 2] <- read.table(file = "Alpha_dependence_1_N_1e+05_FOCuS0_gauss_it_100_.txt", row.names = 1)[,2]
table_alpha [(length(alpha_)+1) : (2*length(alpha_)), 2] <- read.table(file = "Alpha_dependence_2_N_1e+05_FOCuS0_gauss_it_100_.txt", row.names = 1)[,2]
table_alpha [(2*length(alpha_)+1) : (3*length(alpha_)), 2] <- read.table(file = "Alpha_dependence_3_N_1e+05_FOCuS0_gauss_it_100_.txt", row.names = 1)[,2]

#plot---------------------------------------------------------------------|
Plot <- ggplot(table_alpha, aes(x = table_alpha[, 1], y = table_alpha[, 2], color = Dimension)) +
  geom_point(shape=1) +
  geom_line( linetype = "dashed") +
  ylab("Run time, seconds") +
  xlab(expression(x= "Parameter "*alpha* "")) +
  scale_x_continuous(breaks=c(1.0625,1.5,2,3,4,8,16,32,64,128),  trans = "log10",
                     labels = c("1.0625","1.5","2","3","4","8","16","32","64","128")) +
  scale_y_continuous(breaks=c(0,1,5,10,20, 30, 60, 120, 180),  trans = "log10",
                     labels = c('0','1','5','10','20', '30', '60', '120', '180')) +
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text=element_text(size = 13),
        legend.position = "right") +
  scale_color_brewer(palette = 'Paired')

print(Plot)
#save---------------------------------------------------------------------|
pdf(file = "figures/Figure_Gaussian_Model_Run_time_as_the function_of_alpha.pdf",  width = 8, height = 3)
print(Plot)
dev.off()

#ggsave(
#  "Figure_Gaussian_Model_Run_time_as_the function_of_alpha.pdf",
#  plot = Plot,
#  width = 8,
#  height = 3
#)


################################################################################
##   Figure Runtime as the function of p. Gaussian Model                      ##                         
################################################################################
#parameters---------------------------------------------------------------|
Dim <- c (1, 2, 3, 4, 5, 6)
nameDim <- c("MdFOCuS0 (p = 1)",
             "MdFOCuS0 (p = 2)",
             "MdFOCuS0 (p = 3)",
             "MdFOCuS0 (p = 4)",
             "MdFOCuS0 (p = 5)",
             "MdFOCuS0 (p = 6)",
             "MdFOCuS (p = 1)",
             "MdFOCuS (p = 2)",
             "MdFOCuS (p = 3)",
             "MdFOCuS (p = 4)",
             "MdFOCuS (p = 5)",
             "MdFOCuS (p = 6)")

Methods <- c("FOCuS0", "FOCuS")

cost  <- "gauss"

Length <-2^(10:23)
NbSimus <- 100
#read results-------------------------------------------------------------|
NbMethods <- length(Methods)
NbDim <- length(Dim)
NbCosts <- length(cost)
read_file <- NULL
test_res<- data.frame(n = Length, matrix(0, length(Length), 2*NbDim))
colnames(test_res) <- c('n', nameDim)
for (j in  1 : NbMethods) {
  for (k  in 1 : NbDim) {
    read_file <- paste ('Runtime_alpha_2_beta_1_P',
                        Dim[k],
                        Methods[j],
                        cost,
                        "it",
                        NbSimus,
                        '.txt',
                        sep = '_')
    if (Methods[j] == "FOCuS0")
      test_res[, k+1] <- read.table(file = read_file, row.names = 1)[,2]
    else
      test_res[, NbDim+k+1] <- read.table(file = read_file, row.names = 1)[,2]
  }
}

Plot2 <-  ggplot(test_res, aes(n))+
  geom_point(aes(y = test_res[,2],color = "p = 1"), shape = 1, size=2) +
  geom_point(aes(y = test_res[,3],color = "p = 2"), shape = 1, size=2) +
  geom_point(aes(y = test_res[,4],color = "p = 3"), shape = 1, size=2) +
  geom_point(aes(y = test_res[,5],color = "p = 4"), shape = 1, size=2) +
  geom_point(aes(y = test_res[,6],color = "p = 5"), shape = 1, size=2) +
  geom_point(aes(y = test_res[,7],color = "p = 6"), shape = 1, size=2) +
  ylab("Run time, seconds") +
  xlab("Number of data points of time series, n") +
  geom_hline(yintercept = 1200, linetype = "dotted") +
  scale_x_continuous(breaks=c(2^10,2^11,2^12,2^13, 2^14, 2^15, 2^16, 2^17,2^18,2^19,2^20,2^21,2^22,2^23),  trans = "log10",
                     labels = scales::math_format(2^.x, format = log2))+
  scale_y_continuous(breaks=c(1,5,10, 30, 60, 120,300,600,1200),  trans = "log10",
                     labels = c('1','5', '10', '30', '60', '120', "300","600","1200")) +
  theme_bw()+
  scale_color_brewer(palette = 'Paired') +
  geom_line(aes(y = test_res[,2], linetype = "\U03B7\U2081 is known",color = "p = 1")) +
  geom_line(aes(y = test_res[,3], linetype = "\U03B7\U2081 is known",color= "p = 2")) +
  geom_line(aes(y = test_res[,4], linetype = "\U03B7\U2081 is known",color= "p = 3")) +
  geom_line(aes(y = test_res[,5], linetype = "\U03B7\U2081 is known",color= "p = 4")) +
  geom_line(aes(y = test_res[,6], linetype = "\U03B7\U2081 is known",color= "p = 5")) +
  geom_line(aes(y = test_res[,7], linetype = "\U03B7\U2081 is known",color= "p = 6")) +
  geom_line(aes(y = test_res[,8], linetype = "\U03B7\U2081 is unknown", color = "p = 1")) +
  geom_line(aes(y = test_res[,9], linetype = "\U03B7\U2081 is unknown", color= "p = 2")) +
  geom_line(aes(y = test_res[,10], linetype = "\U03B7\U2081 is unknown", color= "p = 3")) +
  geom_line(aes(y = test_res[,11], linetype = "\U03B7\U2081 is unknown", color= "p = 4")) +
  geom_line(aes(y = test_res[,12], linetype = "\U03B7\U2081 is unknown", color= "p = 5")) +
  geom_line(aes(y = test_res[,13], linetype = "\U03B7\U2081 is unknown", color= "p = 6")) +
  labs(colour = "Dimension",  linetype = "MdFOCuS") + 
  ggtitle("Run time (in log-scale)")+
  theme(text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text=element_text(size = 13),
        legend.position = "right")
#SLOPE-------------------------------------------------------------------|
logtest_res <- log10(test_res)
#regression---------------------------------------------------------------|
regres<-list()
alpha<-list()
for (i in 1:10){
  regres[[i]]<- lm(logtest_res[, i+1] ~ logtest_res[, 1])
  alpha[[i]]<-regres[[i]]$coefficients[2]
}

alpha <- unlist(alpha)
alpha <- as.vector(alpha,'numeric')
# P=       1        2       3       4         5       6     
#     1.004549 1.115585 1.179658 1.196108 1.143118 1.201544 
#     1.006225 1.118467 1.182693 1.201839 1.143913 1.197309
#save---------------------------------------------------------------------|
pdf(file = "figures/Figure_Gaussian_Model_Runtime.pdf",  width = 6, height = 4)
print(Plot2)
dev.off()

################################################################################
##   Figure Runtime as the function of p. Poisson Model                      ##                         
################################################################################

#parameters---------------------------------------------------------------|
Dim <- c (1, 2, 3, 4, 5, 6)
nameDim <- c("MdFOCuS0 (p = 1)",
             "MdFOCuS0 (p = 2)",
             "MdFOCuS0 (p = 3)",
             "MdFOCuS0 (p = 4)",
             "MdFOCuS0 (p = 5)",
             "MdFOCuS0 (p = 6)",
             "MdFOCuS (p = 1)",
             "MdFOCuS (p = 2)",
             "MdFOCuS (p = 3)",
             "MdFOCuS (p = 4)",
             "MdFOCuS (p = 5)",
             "MdFOCuS (p = 6)")

Methods <- c("FOCuS0", "FOCuS")

cost  <- "poisson"

Length <-2^(10:23)
NbSimus <- 100

#read results-------------------------------------------------------------|
NbMethods <- length(Methods)
NbDim <- length(Dim)
NbCosts <- length(cost)
read_file <- NULL
test_res<- data.frame(n = Length, matrix(0, length(Length), 2*NbDim))
colnames(test_res) <- c('n', nameDim)
for (j in  1 : NbMethods) {
  for (k  in 1 : NbDim) {
    read_file <- paste ('Runtime_alpha_2_beta_1_P',
                        Dim[k],
                        Methods[j],
                        cost,
                        "it",
                        NbSimus,
                        '.txt',
                        sep = '_')
    if (Methods[j] == "FOCuS0")
      test_res[, k+1] <- read.table(file = read_file, row.names = 1)[,2]
    else
      test_res[, NbDim+k+1] <- read.table(file = read_file, row.names = 1)[,2]
  }
}

Plot3 <-  ggplot(test_res, aes(n)) +
  geom_point(aes(y = test_res[,2],color = "p = 1"), shape = 1, size=2) +
  geom_point(aes(y = test_res[,3],color = "p = 2"), shape = 1, size=2) +
  geom_point(aes(y = test_res[,4],color = "p = 3"), shape = 1, size=2) +
  geom_point(aes(y = test_res[,5],color = "p = 4"), shape = 1, size=2) +
  geom_point(aes(y = test_res[,6],color = "p = 5"), shape = 1, size=2) +
  geom_point(aes(y = test_res[,7],color = "p = 6"), shape = 1, size=2) +
  ylab("Run time, seconds") +
  xlab("Number of data points of time series, n") +
  geom_hline(yintercept = 1200, linetype = "dotted") +
  scale_x_continuous(breaks=c(2^10,2^11,2^12,2^13, 2^14, 2^15, 2^16, 2^17,2^18,2^19,2^20,2^21,2^22,2^23),  trans = "log10",
                     labels = scales::math_format(2^.x, format = log2))+
  scale_y_continuous(breaks=c(1,5,10, 30, 60, 120,300,600, 1200),  trans = "log10",
                     labels = c('1','5', '10', '30', '60', '120', "300","600", "1200")) +
  theme_bw()+
  scale_color_brewer(palette = 'Paired') +
  geom_line(aes(y = test_res[,2], linetype = "\U03B7\U2081 is known",color = "p = 1")) +
  geom_line(aes(y = test_res[,3], linetype = "\U03B7\U2081 is known",color= "p = 2")) +
  geom_line(aes(y = test_res[,4], linetype = "\U03B7\U2081 is known",color= "p = 3")) +
  geom_line(aes(y = test_res[,5], linetype = "\U03B7\U2081 is known",color= "p = 4")) +
  geom_line(aes(y = test_res[,6], linetype = "\U03B7\U2081 is known",color= "p = 5")) +
  geom_line(aes(y = test_res[,7], linetype = "\U03B7\U2081 is known",color= "p = 6")) +
  geom_line(aes(y = test_res[,8], linetype = "\U03B7\U2081 is unknown", color = "p = 1")) +
  geom_line(aes(y = test_res[,9], linetype = "\U03B7\U2081 is unknown",color= "p = 2")) +
  geom_line(aes(y = test_res[,10], linetype = "\U03B7\U2081 is unknown",color= "p = 3")) +
  geom_line(aes(y = test_res[,11], linetype = "\U03B7\U2081 is unknown",color= "p = 4")) +
  geom_line(aes(y = test_res[,12], linetype = "\U03B7\U2081 is unknown",color= "p = 5")) +
  geom_line(aes(y = test_res[,13], linetype = "\U03B7\U2081 is unknown",color= "p = 6")) +
  labs(colour = "Dimension",  linetype = "MdFOCuS") +
  ggtitle("Run time (in log-scale)")+
  theme(text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text=element_text(size = 13),
        legend.position = "right")

#save---------------------------------------------------------------------|
pdf(file = "figures/Figure_Poisson_Model_Runtime.pdf",  width = 6, height = 4)
print(Plot3)
dev.off()


#SLOPE--------------------------------------------------------------------|
logtest_res <- log10(test_res)
#regression---------------------------------------------------------------|
regres<-list()
alpha<-list()
for (i in 1:12){
  regres[[i]]<- lm(logtest_res[, i+1] ~ logtest_res[, 1])
  alpha[[i]]<-regres[[i]]$coefficients[2]
}

alpha <- unlist(alpha)
alpha <- as.vector(alpha,'numeric')
alpha 
# P=       1        2       3       4         5       6     
#    1.013261 1.126261 1.192440 1.210478 1.172368 1.249825
#    1.017782 1.138046 1.212610 1.239983 1.190441 1.253617
################################################################################
########################### END ################################################
################################################################################
