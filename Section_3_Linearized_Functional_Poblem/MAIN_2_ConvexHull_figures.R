################################################################################
##     Figures Convex hull and Random Walk (1D and 2D data)                   ##
################################################################################
#packages-----------------------------------------------------------------|
library(ggplot2)
library(plotly)
library(RColorBrewer)
library(geometry)
#functions----------------------------------------------------------------|
Plot_ConvexHull_2D <- function(data_1D) {
  #colors 
  Col = c("#A6CEE3","#1F78B4","#E31A1C")
  data <- data.frame(t = 1:length(data_1D), St = cumsum(data_1D))
  
  hpts <- chull(data)
  hpts <-c(hpts, hpts[1])
  data_hull <- data[hpts, ]
  plt <- ggplot(data, aes(x = t, y = St)) +
    geom_line(color = Col[2]) +
    geom_path(data = data_hull, aes(x = t, y = St), color = Col[1], linetype = "dashed") +
    geom_point(data = data_hull, aes(x = t, y = St), color = Col[3], size = 1.5) +
    labs(x = "i", y = "S(i)") +
    theme_bw()
  return(plt)
}
#------------------------------------------------------------------------------
Plot_ConvexHull_3D <- function(data_2D) {
  #colors
  Col <- c("#A6CEE3", "#1F78B4", "#E31A1C")
  data <- data.frame(t = 1:nrow(data_2D), 
                     S1t = cumsum(data_2D[, 1]), 
                     S2t = cumsum(data_2D[, 2]))
  #quickhull
  hull <- geometry::convhulln(data)
  #vertices
  pnts <- sort(unique(as.vector(hull)))
  markers <- plot_ly(data[pnts, ], 
                     x = ~t, 
                     y = ~S1t, 
                     z = ~S2t, 
                     type = "scatter3d", mode = "markers", marker = list(size = 3, 
                                                                         color = Col[3]),
                     showlegend = FALSE)
  #RW
  lines <- plot_ly(data, x = ~t, 
                   y = ~S1t, 
                   z = ~S2t, 
                   type = "scatter3d", mode = "lines", 
                   line = list(color = Col[2]),
                   showlegend = FALSE)
  #facets
  for (i in 1:nrow(hull)) {
    lines <- lines %>% add_trace(data = data[hull[i, c(1, 2, 3, 1)], ], 
                                 type = "scatter3d", mode = "lines",
                                 line = list(dash = "dash", color = Col[1]),
                                 showlegend = FALSE)
  }
  #plot
  plot <- subplot(markers, lines) %>% 
    layout(scene = list(
      xaxis = list(title = "i"),
      yaxis = list(title = "S\U00B9(i)"),
      zaxis = list(title = "S\U00B2(i)"), aspectratio = list(x = 1, y = 1, z = 1)
    ))
  return(plot)
}


#----------------------------------------
set.seed(115)
n <- 50
x <- rnorm(n, 0, 1)
plot_2D <- Plot_ConvexHull_2D(x)
pdf(file = "figure/Figure_CHwithRW_1D.pdf",  width = 4, height = 4)
print(plot_2D)
dev.off()
#-------------------------------------

set.seed(17)
y <- rnorm(n, 0, 1)
#---------------------------------------------------------------------------
data_2D <-  cbind (x,y)
plot_3D <- Plot_ConvexHull_3D(data_2D)
print(plot_3D)
#Attention: save in "Figure_CHwithRW_2D.tiff"
################################################################################
##                                     end                                    ##
################################################################################
