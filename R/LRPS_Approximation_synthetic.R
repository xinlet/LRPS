library(Metrics)
library(matrixStats)
library(VARshrink)
#======================================
# Two functions are situations in the synthetic data testing situations
# Where the K values are pre-set.

K = 5 # or 10
#======================================
# Low-Rank Pre-smoothed Data
LRPS_datagenerate_synthetic <- function(y_train) {
  pc <- prcomp(y_train,
               center = TRUE,
               scale. = TRUE)
  y_train <- t(t(pc$x[,1:K] %*% t(pc$rotation[,1:K])) * pc$scale + pc$center)
  y_train <- as.data.frame(y_train) 
}

# LRPS Approximation
LRPS_Approx_synthetic <- function(y_train) {
  pc <- prcomp(y_train,
               center = TRUE,
               scale. = TRUE)
  y_train <- t(t(pc$x[,1:K] %*% t(pc$rotation[,1:K])) * pc$scale + pc$center)
  y_train <- as.data.frame(y_train) 
  sparsevar <- VARshrink(y_train, p=1, type = "none", method = "ridge")
}
