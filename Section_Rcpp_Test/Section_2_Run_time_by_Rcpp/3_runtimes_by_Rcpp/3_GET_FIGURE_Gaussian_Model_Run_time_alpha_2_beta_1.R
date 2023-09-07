################################################################################
##     Figure X Number of Facets and Vertices  of Convex Hull                 ##
##     Empirical results and theoretical  estimations by Stirling numbers     ##                               ##
##                      Gaussian Model (attention : 30 mn)                    ##
################################################################################

#packages-----------------------------------------------------------------|
library(ggplot2)
library(stats)
library(RColorBrewer)
library(ggpubr)
library(scales)

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

Plot <-  ggplot(test_res, aes(n))+
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
  geom_line(aes(y = test_res[,2],    linetype = "MdFOCuS0",color = "p = 1")) +
  geom_line(aes(y = test_res[,3],    linetype = "MdFOCuS0",color= "p = 2")) +
  geom_line(aes(y = test_res[,4],    linetype = "MdFOCuS0",color= "p = 3")) +
  geom_line(aes(y = test_res[,5],    linetype = "MdFOCuS0",color= "p = 4")) +
  geom_line(aes(y = test_res[,6],    linetype = "MdFOCuS0",color= "p = 5")) +
  geom_line(aes(y = test_res[,7],linetype = "MdFOCuS", color = "p = 1")) +
  geom_line(aes(y = test_res[,8],linetype = "MdFOCuS",color= "p = 2")) +
  geom_line(aes(y = test_res[,9],linetype = "MdFOCuS",color= "p = 3")) +
  geom_line(aes(y = test_res[,10],linetype = "MdFOCuS",color= "p = 4")) +
  geom_line(aes(y = test_res[,11], linetype = "MdFOCuS",color= "p = 5")) +
  labs(colour = "Dimension",  linetype = "Method") 

#save---------------------------------------------------------------------|
pdf(file = "Figure_Gaussian_Model_Runtime.pdf",  width = 8, height = 4.5)
print(Plot)
dev.off()


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

#[1] 1.010670 1.142203 1.196649 1.200193 1.186469 1.009385 1.113014 1.206233 1.185324 1.193623
#save---------------------------------------------------------------------|
pdf(file = "Figure_Gaussian_Model_Runtime.pdf",  width = 8, height = 4.5)
print(Plot)
dev.off()
################################################################################
########################### END ################################################
################################################################################

