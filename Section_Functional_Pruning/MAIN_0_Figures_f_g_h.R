################################################################################
##        Linearization of f. functions g and h. FOCuS                        ##
################################################################################


#packages-----------------------------------------------------------------|
library(ggplot2)
library(RColorBrewer)
library(scales)
library(ggpubr)
library(plot3D)
library(cowplot)
library(latex2exp)
library(gridExtra)
library(plotly)
#parameters--------------------------------------------------------------
set.seed(22)
indices <- c(1,2, 3)
#indices <- c(1,5, 6)
n <- 10
x = rnorm(n, 0, 1)
a =  indices 
b =  cumsum(x)
b = b[indices]
c=rep(0,3)

#colors 
Col = c("#A6CEE3", "#B2DF8A", "#FB9A99" )


#functions---------------------------------------------------------------

# A function in Gaussian case (change in mean)
fun_A <- function(mu)  {mu^2}
#f functions
fun_f1 <-  function(mu) { return(a[1] * fun_A(mu) - 2 * b[1] * mu + c[1]) }
fun_f2 <-  function(mu) { a[2] * fun_A(mu) - 2 * b[2] * mu + c[2] }
fun_f3 <-  function(mu) { a[3] * fun_A(mu) - 2 * b[3] * mu + c[3] }
fun_f0max <- function(mu) { pmax(fun_f1(mu),
                                fun_f2(mu)
) }
fun_fmax <- function(mu) { pmax(fun_f1(mu),
                                fun_f2(mu),
                                fun_f3(mu)
) }
#------------------------
#projection functions for f
fun_pr_f1 <- function (mu, nb_f = 3) {
  res <- NULL
  if (nb_f == 2) {res <- as.numeric(!(fun_f1(mu) == fun_f0max(mu)))}
  if (nb_f == 3) {res <- as.numeric(!(fun_f1(mu) == fun_fmax(mu)))} 
  res[ res  == 1] <- NA  
  return (res)
}
fun_pr_f2 <- function (mu, nb_f = 3) {
  if (nb_f == 2) {res <- as.numeric(!(fun_f2(mu) == fun_f0max(mu)))}
  if (nb_f == 3) {res <- as.numeric(!(fun_f2(mu) == fun_fmax(mu)))} 
  res[ res  == 1] <- NA  
  return (res)
}
fun_pr_f3 <- function (mu, nb_f = 3) {
  if (nb_f == 2) {res <- as.numeric(!(fun_f3(mu) == fun_f0max(mu)))}
  if (nb_f == 3) {res <- as.numeric(!(fun_f3(mu) == fun_fmax(mu)))} 
  res[ res  == 1] <- NA  
  return (res)
}
#------------------------------
#Linearization : g functions
fun_g1 <-  function(mu, lambda) { a[1] * lambda - 2 * b[1] * mu + c[1] }
fun_g2 <-  function(mu, lambda) { a[2] * lambda - 2 * b[2] * mu + c[2] }
fun_g3 <-  function(mu, lambda) { a[3] * lambda - 2 * b[3] * mu + c[3] }
fun_gmax <- function(mu, lambda) { pmax (fun_g1(mu,lambda),
                                         fun_g2(mu,lambda),
                                         fun_g3(mu,lambda)
)}
fun_g0max <- function(mu, lambda) { pmax (fun_g1(mu,lambda),
                                         fun_g2(mu,lambda)
)}

#------------------------------
##projection functions for g
fun_pr0_g1 <- function (mu, lambda) {
  res <- as.numeric((fun_g1(mu,lambda) < fun_g0max(mu,lambda)))
  res[ res  == 1] <- NA 
  res[ res  == 0] <- "1"  
  return (res)
}
fun_pr0_g2 <- function (mu, lambda) {
  res <- as.numeric((fun_g2(mu,lambda) < fun_g0max(mu,lambda)))
  res[ res  == 1] <- NA
  res[ res  == 0] <- "2"  
  return (res)
}

fun_pr_g1 <- function (mu, lambda, nb_g = 3) {
  res <- NULL
  if (nb_g == 2) {res <- as.numeric((fun_g1(mu,lambda) < fun_g0max(mu,lambda)))}
  if (nb_g == 3) {res <- as.numeric((fun_g1(mu,lambda) < fun_gmax(mu,lambda)))}
  res[ res  == 1] <- NA 
  res[ res  == 0] <- "1"  
  return (res)
}
fun_pr_g2 <- function (mu, lambda, nb_g = 3) {
  res <- NULL
  if (nb_g == 2) {res <- as.numeric((fun_g2(mu,lambda) < fun_g0max(mu,lambda)))}
  if (nb_g == 3) {res <- as.numeric((fun_g2(mu,lambda) < fun_gmax(mu,lambda)))}
  res[ res  == 1] <- NA
  res[ res  == 0] <- "2"  
  return (res)
}
fun_pr_g3 <- function (mu, lambda) {
  res <- as.numeric((fun_g3(mu,lambda) < fun_gmax(mu,lambda)))
  res[ res  == 1] <- NA  
  res[ res  == 0] <- "3"  
  return (res)
}

####PLOTS---------------------------------
grid_nb <- 250 #Change
mu <- seq(-3,3,length = grid_nb) #change
#table with calculations for f functions

table_f0 <- data.frame(mu = rep(mu,2),fi=c(fun_f1(mu),fun_f2(mu)),
                      i = c(rep("1",grid_nb),rep("2",grid_nb)
                      ))

table_pr_f0max <- data.frame(mu = rep(mu,2),
                            fi=c(fun_pr_f1(mu,2),fun_pr_f2(mu,2)),
                            i = c(rep("1",grid_nb), rep("2",grid_nb)
                            ))



table_f <- data.frame(mu = rep(mu,3),
                      fi=c(fun_f1(mu), fun_f2(mu),fun_f3(mu)),
                      i = c(rep("1",grid_nb), rep("2",grid_nb),rep("3",grid_nb)
                      ))

table_pr_fmax <- data.frame(mu = rep(mu,3),
                            fi=c(fun_pr_f1(mu),fun_pr_f2(mu),fun_pr_f3(mu)),
                            i = c(rep("1",grid_nb), rep("2",grid_nb), rep("3",grid_nb)))


Plot_f0 <- ggplot(table_f0, aes(x = table_f0[, 1], 
                              y = table_f0[, 2], 
                              color = table_f0[, 3])) +
  geom_line(size = 1) +
  ggtitle ("Iteration n = 3") +
  ylab(expression(f[i](mu)))+
  xlab(expression(x = mu)) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values = Col) +
  labs(colour = "index i")

Plot_f <- ggplot(table_f, aes(x = table_f[, 1], 
                              y = table_f[, 2], 
                              color = table_f[, 3])) +
  geom_line(size = 1) +
  ggtitle ("Iteration n = 4") +
  ylab(expression(f[i](mu)))+
  xlab(expression(x = mu)) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values = Col) +
  labs(colour = "index i")

Plot_pr_f0max <- ggplot(table_pr_f0max, aes(x = table_pr_f0max[, 1], 
                                          y = table_pr_f0max[, 2], 
                                          color = table_pr_f0max[, 3])) +
  geom_line(size = 1.25) +
  ylab(expression(y = "")) +
  xlab(expression(x = mu)) +
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text=element_text(size = 13),
        legend.position = "bottom") +
  scale_color_manual(values = Col) + 
  labs(colour = "index i") +
  scale_y_discrete(position = "left", breaks=c(0))


Plot_pr_fmax <- ggplot(table_pr_fmax, aes(x = table_pr_fmax[, 1], 
                                          y = table_pr_fmax[, 2], 
                                          color = table_pr_fmax[, 3])) +
  geom_line(size = 1.25) +
  ylab(expression(y = "")) +
  xlab(expression(x = mu)) +
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text=element_text(size = 13),
        legend.position = "bottom") +
  scale_color_manual(values = Col) + 
  labs(colour = "index i") +
  scale_y_discrete(position = "left", breaks=c(0))

plot_grid(Plot_f0, Plot_pr_f0max, ncol = 1, align = "v")
plot_grid(Plot_f, Plot_pr_fmax, ncol = 1, align = "v")

#save
pdf(file = "Figure_functions_fi_and_fmax_n3.pdf",  width = 3, height = 5)
plot_grid(Plot_f0, Plot_pr_f0max, ncol = 1, align = "v")
dev.off()
#save
pdf(file = "Figure_functions_fi_and_fmax_n4.pdf",  width = 3, height = 5)
plot_grid(Plot_f, Plot_pr_fmax, ncol = 1, align = "v")
dev.off()
#g plots-----------------------------------
#parameter lambda
lambda <- seq(min(fun_A(mu)), max(fun_A(mu)),length = grid_nb)
#calculation of g functions
g1 <- outer(mu, lambda, fun_g1)
g2 <- outer(mu, lambda, fun_g2)
g3 <- outer(mu, lambda, fun_g3)
gmax <- outer(mu, lambda, fun_gmax)

pr_g1 <- outer(mu, lambda, fun_pr_g1)
pr_g2 <- outer(mu, lambda, fun_pr_g2)
pr_g3 <- outer(mu, lambda, fun_pr_g3)

pr0_g1 <- outer(mu, lambda, fun_pr0_g1)
pr0_g2 <- outer(mu, lambda, fun_pr0_g2)


res0_i <- rep(0, grid_nb^2) 
res0_i[!is.na(as.vector(pr0_g1))] <- as.vector(pr0_g1)[!is.na(as.vector(pr0_g1))]
res0_i[!is.na(as.vector(pr0_g2))] <- as.vector(pr0_g2)[!is.na(as.vector(pr0_g2))]

res_i <- rep(0, grid_nb^2) 
res_i[!is.na(as.vector(pr_g1))] <- as.vector(pr_g1)[!is.na(as.vector(pr_g1))]
res_i[!is.na(as.vector(pr_g2))] <- as.vector(pr_g2)[!is.na(as.vector(pr_g2))]
res_i[!is.na(as.vector(pr_g3))] <- as.vector(pr_g3)[!is.na(as.vector(pr_g3))]


table_pr0_gmax <- data.frame(mu = rep(mu, each = grid_nb),
                            lambda =rep(lambda, grid_nb),
                            i = res0_i
)


table_pr_gmax <- data.frame(mu = rep(mu, each = grid_nb),
                            lambda =rep(lambda, grid_nb),
                            i = res_i
)

Col_proj = Col[sort(as.numeric(unique(res_i)))]
Col0_proj = Col[sort(as.numeric(unique(res0_i)))]

#plot: 2D-projection max gi
Plot_pr_gmax <- ggplot(table_pr_gmax, aes(x = table_pr_gmax[, 1], 
                                          y = table_pr_gmax[, 2], 
                                          color = table_pr_gmax[, 3])) +
  geom_point(size = 1, shape = 19) +
  ylab(expression(y =  ""*lambda*"")) +
  xlab(expression(x = ""*mu*"")) +
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text=element_text(size = 13),
        legend.position = "bottom") +
  scale_color_manual(values = Col_proj) + 
  labs(colour = "index i") 
pdf(file = "Figure_function_proj_gmax_n4.pdf",  width = 3, height = 2.5)
print(Plot_pr_gmax)
dev.off()

Plot_pr0_gmax <- ggplot(table_pr0_gmax, aes(x = table_pr0_gmax[, 1], 
                                          y = table_pr0_gmax[, 2], 
                                          color = table_pr0_gmax[, 3])) +
  geom_point(size = 1, shape = 19) +
  ylab(expression(y =  ""*lambda*"")) +
  xlab(expression(x = ""*mu*"")) +
  theme_bw() +
  theme(text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text=element_text(size = 13),
        legend.position = "bottom") +
  scale_color_manual(values = Col0_proj) + 
  labs(colour = "index i") 
pdf(file = "Figure_function_proj_gmax_n3.pdf",  width = 3, height = 2.5)
print(Plot_pr0_gmax)
dev.off()




#plot: gi functions-------------------------------------------
plot_g <- plot_ly() %>% 
  add_trace(
    type = "surface",
    x = mu, 
    y = lambda, 
    z = matrix(g1, nrow = grid_nb, byrow = FALSE),
    colorscale = list(c(0, Col[1]), c(1, Col[1])),
    name = "Surface 1", 
    opacity = 1,
    legendgroup = "Group 1",
    colorbar = list(title = "i = 1")) %>%  
  add_trace(
    type = "surface",
    x = mu, 
    y = lambda, 
    z = matrix(g2, nrow = grid_nb, byrow = FALSE),
    colorscale = list(c(0, Col[2]), c(1, Col[2])),
    name = "Surface 2",  
    opacity = 1,
    legendgroup = "Group 2",  
    colorbar = list(title = "i = 2")  
  )%>%  
  add_trace(
    type = "surface",
    x = mu, 
    y = lambda, 
    z = matrix(g3, nrow = grid_nb, byrow = FALSE),
    colorscale = list(c(0, Col[3]), c(1, Col[3])),
    name = "Surface 3",  
    opacity = 1,
    legendgroup = "Group 3",  
    colorbar = list(title = "i = 3")  
  ) %>%
  layout(scene = list(
    xaxis = list(title = "\U03BC"),
    yaxis = list(title = "\U03BB"),
    zaxis = list(title = "g\U1D62(\U03BC,\U03BB)"),
    aspectratio = list(x = 1, y = 1, z = 1)
  )) %>%
  layout(legend = list(orientation = "h"))  


plot_g
#save in "Figure_functions_g_i_n4.tiff"

plot_g0 <- plot_ly() %>% 
  add_trace(
    type = "surface",
    x = mu, 
    y = lambda, 
    z = matrix(g1, nrow = grid_nb, byrow = FALSE),
    colorscale = list(c(0, Col[1]), c(1, Col[1])),
    name = "Surface 1", 
    opacity = 1,
    legendgroup = "Group 1",
    colorbar = list(title = "i = 1")) %>%  
  add_trace(
    type = "surface",
    x = mu, 
    y = lambda, 
    z = matrix(g2, nrow = grid_nb, byrow = FALSE),
    colorscale = list(c(0, Col[2]), c(1, Col[2])),
    name = "Surface 2",  
    opacity = 1,
    legendgroup = "Group 2",  
    colorbar = list(title = "i = 2")  
  ) %>%
  layout(scene = list(
    xaxis = list(title = "\U03BC"),
    yaxis = list(title = "\U03BB"),
    zaxis = list(title = "g\U1D62(\U03BC,\U03BB)"),
    aspectratio = list(x = 1, y = 1, z = 1)
  )) %>%
  layout(legend = list(orientation = "h"))  


plot_g0
#save in "Figure_functions_g_i_n3.tiff"
#plot h----------------------------------------------------------
kappa <- seq(min(g1,g2,g3), max(g1,g2,g3), 
             length.out = 10)# grid_nb )


library(dplyr)
library(tidyr)

grid <- expand.grid( mu = mu, lambda = lambda, kappa = kappa)
grid$h1 <- with(grid, a[1] * lambda + 2 * b[1] * mu + c[1] - kappa)
grid$h2 <- with(grid, a[2] * lambda + 2 * b[2] * mu + c[2] - kappa)
grid$h3 <- with(grid, a[3] * lambda + 2 * b[3] * mu + c[3] - kappa)

h1  <- grid %>% filter(h1 <= 0)
h2  <- grid %>% filter(h2 <= 0)
h3  <- grid %>% filter(h3 <= 0)


plot1 <- plot_ly(data = h1, 
                 x = ~mu, 
                 y = ~lambda, 
                 z = ~kappa, 
                 type = "scatter3d",
                 mode = "markers", 
                 showlegend = FALSE,
                 marker = list(size = 1, color = Col[1]), opacity = 0.1)
plot2 <- plot_ly(data = h2, 
                 x = ~mu, 
                 y = ~lambda, 
                 z = ~kappa, 
                 type = "scatter3d",
                 mode = "markers", 
                 showlegend = FALSE,
                 marker = list(size = 1, color = Col[2]), opacity = 0.1)
plot3 <- plot_ly(data = h3, 
                 x = ~mu, 
                 y = ~lambda, 
                 z = ~kappa, 
                 type = "scatter3d",
                 mode = "markers", 
                 showlegend = FALSE,
                 marker = list(size = 1, color = Col[3]), opacity = 0.1)


plot_h <- subplot(plot1, plot2, plot3)  %>% 
  layout(scene = list(
    xaxis = list(title = "\U03BC"),
    yaxis = list(title = "\U03BB"),
    zaxis = list(title = "h\U1D62(\U03BA,\U03BC, \U03BB)"), aspectratio = list(x = 1, y = 1, z = 1)
  ))

plot_h
#save in "Figure_functions_h_i_n4.tiff"


#plot h0----------------------------------------------------------



kappa <- seq(min(g1,g2), max(g1,g2), 
             length.out = 10)# grid_nb )

grid0 <- expand.grid( mu = mu, lambda = lambda, kappa = kappa)
grid0$h1 <- with(grid0, a[1] * lambda + 2 * b[1] * mu + c[1] - kappa)
grid0$h2 <- with(grid0, a[2] * lambda + 2 * b[2] * mu + c[2] - kappa)

h01  <- grid0 %>% filter(h1 <= 0)
h02  <- grid0 %>% filter(h2 <= 0)


plot01 <- plot_ly(data = h01, 
                 x = ~mu, 
                 y = ~lambda, 
                 z = ~kappa, 
                 type = "scatter3d",
                 mode = "markers", 
                 showlegend = FALSE,
                 marker = list(size = 1, color = Col[1]), opacity = 0.1)
plot02 <- plot_ly(data = h02, 
                 x = ~mu, 
                 y = ~lambda, 
                 z = ~kappa, 
                 type = "scatter3d",
                 mode = "markers", 
                 showlegend = FALSE,
                 marker = list(size = 1, color = Col[2]), opacity = 0.1)



plot_h0 <- subplot(plot01, plot02)  %>% 
  layout(scene = list(
    xaxis = list(title = "\U03BC"),
    yaxis = list(title = "\U03BB"),
    zaxis = list(title = "h\U1D62(\U03BA,\U03BC, \U03BB)"), aspectratio = list(x = 1, y = 1, z = 1)
  ))

plot_h0
#save in "Figure_functions_h_i_n3.tiff"

