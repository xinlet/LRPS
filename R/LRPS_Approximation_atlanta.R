source("Choose_k.R")
library(Metrics)
library(matrixStats)
library(VARshrink)
#======================================
# The first two functions are situations in the synthetic data testing situations
# Where the K values are pre-set.
# The third and fourth functions are situations where the num of factors are unknown
# So we add an augments K to let user select num of factors
# We can also see the function in Choose_k.R file where it gives a rolling window test on choosing k 
# that K number of factors can explain certain percentage of total variance.
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


# Low-Rank Pre-smoothed Data
LRPS_datagenerate_atlanta <- function(y_train, K) {
  pc <- prcomp(y_train,
               center = TRUE,
               scale. = TRUE)
  y_train <- t(t(pc$x[,1:K] %*% t(pc$rotation[,1:K])) * pc$scale + pc$center)
  y_train <- as.data.frame(y_train) 
}
# LRPS Approximation
LRPS_Approx_atlanta <- function(y_train, K) {
  pc <- prcomp(y_train,
               center = TRUE,
               scale. = TRUE)
  y_train <- t(t(pc$x[,1:K] %*% t(pc$rotation[,1:K])) * pc$scale + pc$center)
  y_train <- as.data.frame(y_train) 
  sparsevar <- VARshrink(y_train, p=1, type = "none", method = "ridge")
}
