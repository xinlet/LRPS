library(Metrics)
library(matrixStats)
library(VARshrink)
#======================================
LRPS_Approx <- function(y_train) {
  pc <- prcomp(y_train,
               center = TRUE,
               scale. = TRUE)
  y_train <- t(t(pc$x[,1:K] %*% t(pc$rotation[,1:K])) * pc$scale + pc$center)
  y_train <- as.data.frame(y_train) 
  
  sparsevar <- VARshrink(y_train, p=1, type = "none", method = "ridge")
}
