################################################################################
##     Figure Empirical Run time as the function of qmin                      ##                             
##                      Gaussian Model (MdFOCuS0, p=2 and 3)                  ##
################################################################################
# Remark------------------------------------------------------------------|
# We use the results obtained in MAIN_XX and MAIN_XXI                     |
#-------------------------------------------------------------------------|

#packages-----------------------------------------------------------------|
library(ggplot2)
library(stats)
library(RColorBrewer)
library(ggpubr)
library(scales)
#parameters---------------------------------------------------------------|
vect_qmin <- 3:12
N <-10^5
NbSimus <- 100
Method <-c("FOCuS0")
Cost <- c("gauss")
#data frame---------------------------------------------------------------|
table_qmin <- data.frame(qmin = c(vect_qmin, vect_qmin, vect_qmin),
                         Runtime = rep(NA, 3*length(vect_qmin)), 
                         Dimension = c( rep("p = 1", length(vect_qmin)),
                                        rep("p = 2", length(vect_qmin)), 
                                        rep("p = 3", length(vect_qmin))
                         ))
colnames(table_qmin) <- c('qmin', 'Run time, seconds', 'Dimension')
table_qmin [1:length(vect_qmin), 2] <- read.table(file = "Qmin_dependence_1_N_1e+05_FOCuS0_gauss_it_100_.txt", row.names = 1)[,2]
table_qmin [(length(vect_qmin)+1) : (2*length(vect_qmin)), 2] <- read.table(file = "Qmin_dependence_2_N_1e+05_FOCuS0_gauss_it_100_.txt", row.names = 1)[,2]
table_qmin [(2*length(vect_qmin)+1) : (3*length(vect_qmin)), 2] <- read.table(file = "Qmin_dependence_3_N_1e+05_FOCuS0_gauss_it_100_.txt", row.names = 1)[,2]

#plot---------------------------------------------------------------------|
Plot <- ggplot(table_qmin, aes(x = table_qmin[, 1], y = table_qmin[, 2], color = Dimension)) +
  geom_point(shape=1) +
  geom_line( linetype = "dashed") +
  ylab("Run time, seconds") +
  xlab("Parameter qmin") +
  scale_x_continuous(breaks=3:12, 
                     labels = c("3","4","5","6","7","8","9","10", "11", "12")) +
  scale_y_continuous(breaks= 0:22, 
                     labels = c('0','1','2','3','4','5','6','7','8','9','10','11', '12', '13', '14', '15','16', '17', '18', '19', '20', '21', '22')) +
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text=element_text(size = 13),
        legend.position = "right") +
  scale_color_brewer(palette = 'Paired') 

print(Plot)
#save---------------------------------------------------------------------|
pdf(file = "figure/Figure_Dyadic_GM_Run_time_qmin_N10000.pdf",  width = 4, height = 4)
print(Plot)
dev.off()

################################################################################
##   Figure Runtime as the function of p. Dyadic FOCuS Gaussian Model         ##                         
################################################################################
#parameters---------------------------------------------------------------|
Dim <- c (1, 2, 3)
qmin_opt <- c(6,7,8)
nameDim <- c("dyadic MdFOCuS0 (p = 1)",
             "dyadic MdFOCuS0 (p = 2)",
             "dyadic MdFOCuS0 (p = 3)",
             "dyadic MdFOCuS (p = 1)",
             "dyadic MdFOCuS (p = 2)",
             "dyadic MdFOCuS (p = 3)",
             "MdFOCuS0 (p = 1)",
             "MdFOCuS0 (p = 2)",
             "MdFOCuS0 (p = 3)",
             "MdFOCuS (p = 1)",
             "MdFOCuS (p = 2)",
             "MdFOCuS (p = 3)"
)

Methods <- c("FOCuS0", "FOCuS")

cost  <- "gauss"

Length <-2^(10:23)
NbSimus <- 100
#read results-------------------------------------------------------------|
NbMethods <- length(Methods)
NbDim <- length(Dim)
NbCosts <- length(cost)
read_file <- NULL
test_res <- data.frame(n = Length, matrix(0, length(Length), 4*NbDim))
colnames(test_res) <- c('n', nameDim)
#read Dyadic MdFOCuS
for (j in  1 : NbMethods) {
  for (k  in 1 : NbDim) {
    read_file <- paste ('Runtime_qmin',qmin_opt[k] ,'P',
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
#read MdFOCuS
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
      test_res[, 2 * NbDim + k + 1] <- read.table(file = read_file, row.names = 1)[,2]
    else
      test_res[, 3 * NbDim + k + 1] <- read.table(file = read_file, row.names = 1)[,2]
  }
}


PlotDyadicvsClassic <-  ggplot(test_res, aes(n))+
  ylab("Run time, seconds") +
  xlab("Number of data points of time series, n") +
  geom_hline(yintercept = 180, linetype = "dotted") +
  scale_x_continuous(breaks=c(2^10,2^11,2^12,2^13, 2^14, 2^15, 2^16, 2^17,2^18,2^19,2^20,2^21,2^22,2^23),  trans = "log10",
                     labels = scales::math_format(2^.x, format = log2))+
  scale_y_continuous(breaks=c(1,5,10,20, 30, 60, 120, 180),  trans = "log10",
                     labels = c('1','5', '10','20', '30', '60', '120', "180")) +
  theme_bw()+
  scale_color_brewer(palette = 'Paired') +
  #Dyadic
  geom_point(aes(y = test_res[,2],color = "p = 1"), shape = 1, size=1) +
  geom_point(aes(y = test_res[,3],color = "p = 2"), shape = 1, size=1) +
  geom_point(aes(y = test_res[,4],color = "p = 3"), shape = 1, size=1) +
  #Classic 
  geom_point(aes(y = test_res[,8],color = "p = 1"), shape = 1, size=1) +
  geom_point(aes(y = test_res[,9],color = "p = 2"), shape = 1, size=1) +
  geom_point(aes(y = test_res[,10],color = "p = 3"), shape = 1, size=1) +
  #Dyadic
  geom_line(aes(y = test_res[,2], linetype = "dyadic MdFOCuS0",color = "p = 1")) +
  geom_line(aes(y = test_res[,3], linetype = "dyadic MdFOCuS0",color= "p = 2")) +
  geom_line(aes(y = test_res[,4], linetype = "dyadic MdFOCuS0",color= "p = 3")) +
#  geom_line(aes(y = test_res[,5], linetype = "dyadic MdFOCuS", color = "p = 1")) +
#  geom_line(aes(y = test_res[,6], linetype = "dyadic MdFOCuS", color= "p = 2")) +
#  geom_line(aes(y = test_res[,7], linetype = "dyadic MdFOCuS", color= "p = 3")) +
  #Classic
  geom_line(aes(y = test_res[,8], linetype = "MdFOCuS0",color = "p = 1")) +
  geom_line(aes(y = test_res[,9], linetype = "MdFOCuS0",color= "p = 2")) +
  geom_line(aes(y = test_res[,10], linetype = "MdFOCuS0",color= "p = 3")) +
#  geom_line(aes(y = test_res[,11], linetype = "MdFOCuS", color = "p = 1")) +
#  geom_line(aes(y = test_res[,12], linetype = "MdFOCuS", color= "p = 2")) +
#  geom_line(aes(y = test_res[,13], linetype = "MdFOCuS", color= "p = 3")) +
  labs(colour = "Dimension",  linetype = "Method") + 
  theme(text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text=element_text(size = 13),
        legend.position = "right")
##save---------------------------------------------------------------------|
pdf(file = "figure/Figure_DyadicvsClassic_GM_Runtime_p.pdf",  width = 6, height = 4)
print(PlotDyadicvsClassic )
dev.off()
#------------------------------------------------------------------------------
ratio_test_res <- test_res
ratio_test_res[,2:ncol(ratio_test_res)] <- log10(ratio_test_res[,2:ncol(ratio_test_res)]/ratio_test_res[,1])
ratio_test_res[,1] <- log10(log10(ratio_test_res[,1]))

#regression log(T/n) = beta log(log(n)) ---------------------------------------|
#Dyadic
regres_dyadic <- list()
beta_dyadic <- list()
for (i in 1:6){
  regres_dyadic[[i]]<- lm(ratio_test_res[, i+1] ~ ratio_test_res[, 1])
  beta_dyadic[[i]]<-regres_dyadic[[i]]$coefficients[2]
}

beta_dyadic <- unlist(beta_dyadic)
beta_dyadic <- as.vector(beta_dyadic,'numeric')
beta_dyadic
#[1]FOCuS0: 1.149658 2.243110 3.203183 FOCuS: 1.176664 2.228784 3.214095

#Classic
regres_classic <- list()
beta_classic <- list()
for (i in 1:6) {
  regres_classic[[i]]<- lm(ratio_test_res[, 6 + i+1] ~ ratio_test_res[, 1])
  beta_classic[[i]]<-regres_classic[[i]]$coefficients[2]
}

beta_classic <- unlist(beta_classic)
beta_classic <- as.vector(beta_classic,'numeric')
beta_classic
#[1]FOCuS0:  0.3788628 1.6531348 2.2241920  FOCuS: 0.4553198 1.6657652 2.2746424


#RATIO
PlotRatioDyadicvsClassic0 <-  ggplot(ratio_test_res, aes(n))+
  ylab("log(Runtime(n)/n)") +
  xlab("log(log(n))") +
  theme_bw()+
  scale_color_brewer(palette = 'Paired') +
  #Dyadic
  geom_point(aes(y = ratio_test_res[,2],color = "p = 1"), shape = 1, size=1) +
  geom_point(aes(y = ratio_test_res[,3],color = "p = 2"), shape = 1, size=1) +
  geom_point(aes(y = ratio_test_res[,4],color = "p = 3"), shape = 1, size=1) +
  #Classic 
  geom_point(aes(y = ratio_test_res[,8],color = "p = 1"), shape = 1, size=1) +
  geom_point(aes(y = ratio_test_res[,9],color = "p = 2"), shape = 1, size=1) +
  geom_point(aes(y = ratio_test_res[,10],color = "p = 3"), shape = 1, size=1) +
  #Dyadic
  geom_line(aes(y = ratio_test_res[,2], linetype = "dyadic MdFOCuS0",color = "p = 1")) +
  geom_line(aes(y = ratio_test_res[,3], linetype = "dyadic MdFOCuS0",color= "p = 2")) +
  geom_line(aes(y = ratio_test_res[,4], linetype = "dyadic MdFOCuS0",color= "p = 3")) +
 # geom_line(aes(y = ratio_test_res[,5], linetype = "dyadic MdFOCuS", color = "p = 1")) +
#  geom_line(aes(y = ratio_test_res[,6], linetype = "dyadic MdFOCuS", color= "p = 2")) +
#  geom_line(aes(y = ratio_test_res[,7], linetype = "dyadic MdFOCuS", color= "p = 3")) +
  #Classic
  geom_line(aes(y = ratio_test_res[,8], linetype = "MdFOCuS0",color = "p = 1")) +
  geom_line(aes(y = ratio_test_res[,9], linetype = "MdFOCuS0",color= "p = 2")) +
  geom_line(aes(y = ratio_test_res[,10], linetype = "MdFOCuS0",color= "p = 3")) +
#  geom_line(aes(y = ratio_test_res[,11], linetype = "MdFOCuS", color = "p = 1")) +
 # geom_line(aes(y = ratio_test_res[,12], linetype = "MdFOCuS", color= "p = 2")) +
 # geom_line(aes(y = ratio_test_res[,13], linetype = "MdFOCuS", color= "p = 3")) +
  labs(colour = "Dimension",  linetype = "Method") + 
  theme(text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text=element_text(size = 13),
        legend.position = "right")
##save---------------------------------------------------------------------|
#pdf(file = "Figure_DyadicvsClassic0_GM_Ratio_Runtime_p.pdf",  width = 6, height = 4)
#print(PlotRatioDyadicvsClassic0)
#dev.off()






