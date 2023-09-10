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
table_alpha <- data.frame(alpha = c(alpha_, alpha_),
                                 Runtime = rep(NA, 2*length(alpha_)), Dimension = c(rep("p = 2", length(alpha_)), rep("p = 3", length(alpha_))))
  colnames(table_alpha) <- c('\u03b1', 'Run time, seconds', 'Dimension')
table_alpha [1:length(alpha_), 2] <- read.table(file = "Alpha_dependence_2_N_1e+05_FOCuS0_gauss_it_100_.txt", row.names = 1)[,2]
table_alpha [(length(alpha_)+1) : (2*length(alpha_)), 2] <- read.table(file = "Alpha_dependence_3_N_1e+05_FOCuS0_gauss_it_100_.txt", row.names = 1)[,2]

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
    scale_color_manual(values = c("#1F78B4","#B2DF8A"))
  
print(Plot)
#save---------------------------------------------------------------------|
pdf(file = "Figure_Gaussian_Model_Run_time_as_the function_of_alpha.pdf",  width = 8, height = 3)
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
Dim <- c (1, 2, 3, 4, 5)
nameDim <- c("MdFOCuS0 (p = 1)",
             "MdFOCuS0 (p = 2)",
             "MdFOCuS0 (p = 3)",
             "MdFOCuS0 (p = 4)",
             "MdFOCuS0 (p = 5)",
             "MdFOCuS (p = 1)",
             "MdFOCuS (p = 2)",
             "MdFOCuS (p = 3)",
             "MdFOCuS (p = 4)",
             "MdFOCuS (p = 5)")

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
  ylab("Run time, seconds") +
  xlab("Number of data points of time series, n") +
  geom_hline(yintercept = 180, linetype = "dotted") +
  scale_x_continuous(breaks=c(2^10,2^11,2^12,2^13, 2^14, 2^15, 2^16, 2^17,2^18,2^19,2^20,2^21,2^22,2^23),  trans = "log10",
                     labels = scales::math_format(2^.x, format = log2))+
  scale_y_continuous(breaks=c(1,5,10,20, 30, 60, 120, 180),  trans = "log10",
                     labels = c('1','5', '10','20', '30', '60', '120', "180")) +
  theme_bw()+
  scale_color_brewer(palette = 'Paired') +
  geom_line(aes(y = test_res[,2], linetype = "MdFOCuS0",color = "p = 1")) +
  geom_line(aes(y = test_res[,3], linetype = "MdFOCuS0",color= "p = 2")) +
  geom_line(aes(y = test_res[,4], linetype = "MdFOCuS0",color= "p = 3")) +
  geom_line(aes(y = test_res[,5], linetype = "MdFOCuS0",color= "p = 4")) +
  geom_line(aes(y = test_res[,6], linetype = "MdFOCuS0",color= "p = 5")) +
  geom_line(aes(y = test_res[,7], linetype = "MdFOCuS", color = "p = 1")) +
  geom_line(aes(y = test_res[,8], linetype = "MdFOCuS", color= "p = 2")) +
  geom_line(aes(y = test_res[,9], linetype = "MdFOCuS", color= "p = 3")) +
  geom_line(aes(y = test_res[,10], linetype = "MdFOCuS", color= "p = 4")) +
  geom_line(aes(y = test_res[,11], linetype = "MdFOCuS", color= "p = 5")) +
  labs(colour = "Dimension",  linetype = "Method") + 
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
#1.038808 1.157248 1.222801 1.245924 1.137845 1.045247 1.158712 1.228089 1.244636 1.140385
#save---------------------------------------------------------------------|
pdf(file = "Figure_Gaussian_Model_Runtime.pdf",  width = 8, height = 4)
print(Plot2)
dev.off()
  
################################################################################
##   Figure Runtime as the function of p. Poisson Model                      ##                         
################################################################################
  
#parameters---------------------------------------------------------------|
Dim <- c (1, 2, 3, 4, 5)
nameDim <- c("MdFOCuS0 (p = 1)",
             "MdFOCuS0 (p = 2)",
             "MdFOCuS0 (p = 3)",
             "MdFOCuS0 (p = 4)",
             "MdFOCuS0 (p = 5)",
             "MdFOCuS (p = 1)",
             "MdFOCuS (p = 2)",
             "MdFOCuS (p = 3)",
             "MdFOCuS (p = 4)",
             "MdFOCuS (p = 5)")

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
  ylab("Run time, seconds") +
  xlab("Number of data points of time series, n") +
  geom_hline(yintercept = 180, linetype = "dotted") +
  scale_x_continuous(breaks=c(2^10,2^11,2^12,2^13, 2^14, 2^15, 2^16, 2^17,2^18,2^19,2^20,2^21,2^22,2^23),  trans = "log10",
                     labels = scales::math_format(2^.x, format = log2))+
  scale_y_continuous(breaks=c(1,5,10,20, 30, 60, 120, 180),  trans = "log10",
                     labels = c('1','5', '10','20', '30', '60', '120', "180")) +
  theme_bw()+
  scale_color_brewer(palette = 'Paired') +
  geom_line(aes(y = test_res[,2], linetype = "MdFOCuS0",color = "p = 1")) +
  geom_line(aes(y = test_res[,3], linetype = "MdFOCuS0",color= "p = 2")) +
  geom_line(aes(y = test_res[,4], linetype = "MdFOCuS0",color= "p = 3")) +
  geom_line(aes(y = test_res[,5], linetype = "MdFOCuS0",color= "p = 4")) +
  geom_line(aes(y = test_res[,6], linetype = "MdFOCuS0",color= "p = 5")) +
  geom_line(aes(y = test_res[,7], linetype = "MdFOCuS", color = "p = 1")) +
  geom_line(aes(y = test_res[,8], linetype = "MdFOCuS",color= "p = 2")) +
  geom_line(aes(y = test_res[,9], linetype = "MdFOCuS",color= "p = 3")) +
  geom_line(aes(y = test_res[,10], linetype = "MdFOCuS",color= "p = 4")) +
  geom_line(aes(y = test_res[,11], linetype = "MdFOCuS",color= "p = 5")) +
  labs(colour = "Dimension",  linetype = "Method") +
  theme(text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text=element_text(size = 13),
        legend.position = "right")

#save---------------------------------------------------------------------|
pdf(file = "Figure_Poisson_Model_Runtime.pdf",  width = 8, height = 4)
print(Plot3)
dev.off()


#SLOPE--------------------------------------------------------------------|
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
alpha 
# [1] 1.044838 1.161654 1.244333 1.259823 1.163840 1.057966 1.182293 1.264088
#[9] 1.286759 1.193310
################################################################################
########################### END ################################################
################################################################################


