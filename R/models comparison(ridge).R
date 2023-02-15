library(vars)
library(tseries)
library(tidyverse)
library(forecast)
library(dplyr)
library(tidyr)
library(sparsevar)
library(Metrics)
library(matrixStats)
library(VARshrink)

#realdata=====================================
setwd("~/UoB/Code/Data_Atlanta fMRI")
filenames <- list.files(path = "~/UoB/Code/Data_Atlanta fMRI",  
                        ("csv"))
K = 10
rollingwindow = 120

#simulation Alex==========================
setwd("~/UoB/Code/Xinle_Models/Pre-smoothing/Simulated DFM new/Alex/k=10/rho=0.7/d=0.8 tau=0.8")
filenames <- list.files(path = "~/UoB/Code/Xinle_Models/Pre-smoothing/Simulated DFM new/Alex/k=10/rho=0.7/d=0.8 tau=0.8",  
                        ("csv"))
K = 10
rollingwindow = 50

#Simplified datasets====
setwd("~/UoB/Code/Xinle_Models/Pre-smoothing/Simulated DFM new/Specific Section/k=5/rho=0.3")
filenames <- list.files(path = "~/UoB/Code/Xinle_Models/Pre-smoothing/Simulated DFM new/Specific Section/k=5/rho=0.3",  
                        ("csv"))
K = 5
rollingwindow = 50
#filenames=============================
filenames <- filenames[1:length(filenames)]

#load====
load("result.RData")

##################################
#1. PCA and Sprase VAR (new generating training data)====================
start.time1 <- Sys.time()
analyzepcasparsevar <- function(filename) {
  y <- read.csv(file = filename, header = TRUE)
  y <- data.matrix(y)
  results = array(0, dim = c(2, (nrow(y) - rollingwindow - 1))) # where to store the results
  for(i in 1:(nrow(y) - rollingwindow - 1)) {
    y_train <- y[i: (i-1+rollingwindow),]   #training data
    y_test <- y[(i+ rollingwindow),]
    
    pc <- prcomp(y_train,
                 center = TRUE,
                 scale. = TRUE)
    y_train <- t(t(pc$x[,1:K] %*% t(pc$rotation[,1:K])) * pc$scale + pc$center)
    y_train <- as.data.frame(y_train) #new train 
    # fit sparseVAR
    sparsevar <- VARshrink(y_train, p=1, type = "none", method = "ridge")
    pred_var <- predict(sparsevar, n.ahead = 1)
    
    pred_var_value <- lapply(pred_var$fcst, `[[`, 1)
    pred_var_value <- as.data.frame(pred_var_value)
    pred_var_num <- as.numeric(unlist(pred_var_value))

    y_test_num <- as.numeric(unlist(y_test))
    rmse_result <- rmse(y_test_num, pred_var_num)
    mae_result <- mae(y_test_num, pred_var_num)
    results[,i] <- c(rmse_result,mae_result)
  }
  Meanresultsind <- c(rowMeans(results))
  return(Meanresultsind)
}

# Build a container matrix
# Container matrix: rows = MSE,RMSE...; cols = Files 
result_pcasparsevar <- matrix(NA, 2, length(filenames), dimnames = list(c( "rmse","mae"),
                                                                        c(as.list(filenames))))

# iterations for each file
for (i in 1:length(filenames)) {
  filenamesuse <- filenames[i]
  output <- analyzepcasparsevar(filenamesuse) 
  result_pcasparsevar[,i] = output
  
}
mean1 <- rowMeans(result_pcasparsevar)
end.time1 <- Sys.time()
time1 <- end.time1 - start.time1
time1

#############################################################
#2. DFM PCA and SVAR =========================

start.time2 <- Sys.time()
analyzepcaSVAR <- function(filename) {
  y <- read.csv(file = filename, header = TRUE)
  y <- data.matrix(y)
  results = array(0, dim = c(2, (nrow(y) - rollingwindow - 1))) # where to store the results
  
  for(i in 1:(nrow(y) - rollingwindow - 1)) {
    y_train <- y[i: (i-1+rollingwindow),]   #training data
    y_test <- y[(i+rollingwindow),]
    
    pc <- prcomp(y_train,
                 center = TRUE,
                 scale. = TRUE)
    
    est <- VARshrink(pc$x[,1:K], p=1, type = "none", method = "ridge")
    pred_var <- predict(est, n.ahead = 1)
    pred_var_value <- lapply(pred_var$fcst, `[[`, 1)
    pred_var_value <- as.data.frame(pred_var_value)
    pred_var_value_num <- as.numeric(unlist(pred_var_value))
    
    y_new <- t(t(rbind(pc$x[,1:K],pred_var_value_num) %*% t(pc$rotation[,1:K])) 
               * pc$scale + pc$center)
    y_new <- as.data.frame(y_new)
    pred_var <- tail(y_new, n = 1)
    
    pred_var_num <- as.numeric(unlist(pred_var))
    y_test_num <- as.numeric(unlist(y_test))
    
    rmse_result <- rmse(y_test_num, pred_var_num)
    mae_result <- mae(y_test_num, pred_var_num)
    results[,i] <- c(rmse_result,mae_result)
  }
  
  Meanresultsind <- c(rowMeans(results))
  return(Meanresultsind)
}


result_pcaSVAR <- matrix(NA, 2, length(filenames), dimnames = list(c("rmse","mae"),
                                                                   c(as.list(filenames))))

# iterations for each file
for (i in 1:length(filenames)) {
  filenamesuse <- filenames[i]
  output <- analyzepcaSVAR(filenamesuse) 
  result_pcaSVAR[,i] = output
  
}
mean2 <- rowMeans(result_pcaSVAR)
end.time2 <- Sys.time()
time2 <- end.time2 - start.time2

#############################################
# 3. DFM PCA + VAR========================
start.time3 <- Sys.time()
analyzepcaVAR <- function(filename) {
  y <- read.csv(file = filename, header = TRUE)
  y <- data.matrix(y)
  results = array(0, dim = c(2, (nrow(y) - rollingwindow - 1))) # where to store the results
  
  for(i in 1:(nrow(y) - rollingwindow - 1)) {
    y_train <- y[i: (i-1+rollingwindow),]   #training data
    y_test <- y[(i+rollingwindow),]
    
    pc <- prcomp(y_train,
                 center = TRUE,
                 scale. = TRUE)
    
    est <- VAR(pc$x[,1:K], p=1, type ="none")
    pred_var <- predict(est, n.ahead = 1)
    pred_var_value <- lapply(pred_var$fcst, `[[`, 1)
    pred_var_value_num <- as.numeric(unlist(pred_var_value))
    
    y_new <- t(t(rbind(pc$x[,1:K],pred_var_value_num) %*% t(pc$rotation[,1:K])) 
               * pc$scale + pc$center)
    y_new <- as.data.frame(y_new)
    pred_var <- tail(y_new, n = 1)
    
    pred_var_num <- as.numeric(unlist(pred_var))
    y_test_num <- as.numeric(unlist(y_test))
    
    rmse_result <- rmse(y_test_num, pred_var_num)
    mae_result <- mae(y_test_num, pred_var_num)
    result = c(rmse_result,mae_result)
    return(result)
  }
  
  Meanresultsind <- c(rowMeans(results))
  return(Meanresultsind)
}

result_pcaVAR <- matrix(NA, 2, length(filenames), dimnames = list(c( "rmse","mae"),
                                                                  c(as.list(filenames))))

for (i in 1:length(filenames)) {
  filenamesuse <- filenames[i]
  output <- analyzepcaVAR(filenamesuse) 
  result_pcaVAR[,i] = output
}
mean3 <- rowMeans(result_pcaVAR)
end.time3 <- Sys.time()
time3 <- end.time3 - start.time3

#####################################
# 4. VAR===========================
start.time4 <- Sys.time()
analyzevar <- function(filename) {
  y <- read.csv(file = filename, header = TRUE)
  y <- data.matrix(y)
  results = array(0, dim = c(2, (nrow(y) - rollingwindow - 1))) # where to store the results
  
  for(i in 1:(nrow(y) - rollingwindow - 1)) {
    y_train <- y[i: (i-1+rollingwindow),]   #training data
    y_test <- y[(i+rollingwindow),]
    
    # fit VAR
    est <- VAR(y_train, p=1, type ="trend") #VAR
    pred_var <- predict(est, n.ahead = 1)
    
    pred_var_value <- lapply(pred_var$fcst, `[[`, 1)
    pred_var_value <- as.data.frame(pred_var_value)
    pred_var_num <- as.numeric(unlist(pred_var_value))
    y_test_num <- as.numeric(unlist(y_test))

    rmse_result <- rmse(y_test_num, pred_var_num)
    mae_result <- mae(y_test_num, pred_var_num)
    result = c(rmse_result,mae_result)
    return(result)
  }
  Meanresultsind <- c(rowMeans(results))
  return(Meanresultsind)
}
# Build a container matrix
# Container matrix: rows = MSE,RMSE...; cols = Files 
result_var <- matrix(NA, 2, length(filenames), dimnames = list(c( "rmse","mae"),
                                                               c(as.list(filenames))))

# iterations for each file
for (i in 1:length(filenames)) {
  filenamesuse <- filenames[i]
  output <- analyzevar(filenamesuse) 
  result_var[,i] = output
  
}
mean4 <- rowMeans(result_var)
end.time4 <- Sys.time()
time4 <- end.time4 - start.time4

##########################################
# 5. low rank + idiosyncratic SVAR===================

start.time5 <- Sys.time()
analyzelowrankSVAR <- function(filename) {
  y <- read.csv(file = filename, header = TRUE)
  y <- data.matrix(y)
  results = array(0, dim = c(2, (nrow(y) - rollingwindow - 1))) # where to store the results
  
  for(i in 1:(nrow(y) - rollingwindow - 1)) {
    y_train <- y[i: (i-1+rollingwindow),]   #training data
    y_test <- y[(i+rollingwindow),]
    
    pc <- prcomp(y_train,
                 center = TRUE,
                 scale. = TRUE)
    y_k <- t(t(pc$x[,1:K] %*% t(pc$rotation[,1:K])) * pc$scale + pc$center)
    y_res <- y_train - y_k
    
    loads <- pc$rotation[,1:K]
    pred_pc <- (loads*y_test)
    
    sparsevar <- VARshrink(y_res, p=1, type = "none", method = "ridge")
    pred_var <- predict(sparsevar, n.ahead = 1)
    pred_var_value <- lapply(pred_var$fcst, `[[`, 1)
    pred_var_value <- as.data.frame(pred_var_value)
    pred_var_value <- pred_var_value + pred_pc
    pred_var_value_num <- as.numeric(unlist(pred_var_value))
    y_test_num <- as.numeric(unlist(y_test))

    rmse_result <- rmse(y_test_num, pred_var_value_num)
    mae_result <- mae(y_test_num, pred_var_value_num)
    results[,i] <- c(rmse_result,mae_result)
  }
  
  Meanresultsind <- c(rowMeans(results))
  return(Meanresultsind)
}


result_lowrankSVAR <- matrix(NA, 2, length(filenames), dimnames = list(c( "rmse","mae"),
                                                                       c(as.list(filenames))))

# iterations for each file
for (i in 1:length(filenames)) {
  filenamesuse <- filenames[i]
  output <- analyzelowrankSVAR(filenamesuse) 
  result_lowrankSVAR[,i] = output
  
}
mean5 <- rowMeans(result_lowrankSVAR)
end.time5 <- Sys.time()
time5 <- end.time5 - start.time5

######################################
# 6. sparse VAR==============================
start.time6 <- Sys.time()
analyzesparsevar <- function(filename) {
  y <- read.csv(file = filename, header = TRUE)
  y <- data.matrix(y)
  results = array(0, dim = c(2, (nrow(y) - rollingwindow - 1))) # where to store the results
  
  for(i in 1:(nrow(y) - rollingwindow - 1)) {
    y_train <- y[i: (i-1+rollingwindow),]   #training data
    y_test <- y[(i+rollingwindow),]
    
    
    # fit sparseVAR
    sparsevar <- VARshrink(y_train, p=1, type = "none", method = "ridge") 
    pred_var <- predict(sparsevar, n.ahead = 1)
    pred_var_value <- lapply(pred_var$fcst, `[[`, 1)
    pred_var_value <- as.data.frame(pred_var_value)
    pred_var_num <- as.numeric(unlist(pred_var_value))

    y_test_num <- as.numeric(unlist(y_test))
    rmse_result <- rmse(y_test_num, pred_var_num)
    mae_result <- mae(y_test_num, pred_var_num)
    result = c(rmse_result,mae_result)
    return(result)
  }
  Meanresultsind <- c(rowMeans(results))
  return(Meanresultsind)
}

# Build a container matrix
# Container matrix: rows = MSE,RMSE...; cols = Files 
result_sparsevar <- matrix(NA, 2, length(filenames), dimnames = list(c( "rmse","mae"),
                                                                     c(as.list(filenames))))

# iterations for each file
for (i in 1:length(filenames)) {
  filenamesuse <- filenames[i]
  output <- analyzesparsevar(filenamesuse) 
  result_sparsevar[,i] = output
  
}
mean6 <- rowMeans(result_sparsevar)
end.time6 <- Sys.time()
time6 <- end.time6 - start.time6

##################################
# 7. sparseVAR lasso====
start.time7 <- Sys.time()
sparsevar_lasso <- function(filename) {
  y <- read.csv(file = filename, header = TRUE)
  y <- data.matrix(y)
  results = array(0, dim = c(2, (nrow(y) - rollingwindow - 1))) # where to store the results
  
  for(i in 1:(nrow(y) - rollingwindow - 1)) {
    y_train <- y[i: (i-1+rollingwindow),]   #training data
    y_test <- y[(i+rollingwindow),]
    
    
    # fit sparseVAR
    sparsevar <- fitVAR(y_train, p=1, penalty = "ENET", method = "cv") #default = lasso
    
    pred_var <- computeForecasts(sparsevar, 1)
    pred_var <- as.data.frame(t(pred_var))
    pred_var <- pred_var[1,]
    
    pred_var_num <- as.numeric(unlist(pred_var))
    y_test_num <- as.numeric(unlist(y_test))
    
    rmse_result <- rmse(y_test_num, pred_var_num)
    mae_result <- mae(y_test_num, pred_var_num)
   
    result = c(rmse_result,mae_result)
    return(result)
  }
  Meanresultsind <- c(rowMeans(results))
  return(Meanresultsind)
}

# Build a container matrix
# Container matrix: rows = MSE,RMSE...; cols = Files 
result_sparsevar_lasso <- matrix(NA, 2, length(filenames), dimnames = list(c("rmse","mae"),
                                                                     c(as.list(filenames))))

# iterations for each file
for (i in 1:length(filenames)) {
  filenamesuse <- filenames[i]
  output <- sparsevar_lasso(filenamesuse) 
  result_sparsevar_lasso[,i] = output
  
}
mean7 <- rowMeans(result_sparsevar_lasso)
end.time7 <- Sys.time()
time7 <- end.time7 - start.time7

#############################################
# 8. DFM sparseVAR lasso====
start.time8 <- Sys.time()
dfmsvar_lasso <- function(filename) {
  y <- read.csv(file = filename, header = TRUE)
  results = array(0, dim = c(2, (nrow(y) - rollingwindow - 1))) # where to store the results
  
  for(i in 1:(nrow(y) - rollingwindow - 1)) {
    y_train <- y[i: (i-1+rollingwindow),]   #training data
    y_test <- y[(i+rollingwindow),]
    
    pc <- prcomp(y_train,
                 center = TRUE,
                 scale. = TRUE)
    
    est <- fitVAR(pc$x[,1:K], p=1, penalty = "ENET", method = "cv")
    pred_var <- computeForecasts(est, 1)
    pred_var_value <- as.data.frame(t(pred_var))
    pred_var_value <- pred_var[1,]
    pred_var_value_num <- as.numeric(unlist(pred_var_value))
    
    y_new <- t(t(rbind(pc$x[,1:K],pred_var_value_num) %*% t(pc$rotation[,1:K])) 
               * pc$scale + pc$center)
    y_new <- as.data.frame(y_new)
    pred_var <- tail(y_new, n = 1)
    
    pred_var_num <- as.numeric(unlist(pred_var))
    y_test_num <- as.numeric(unlist(y_test))
    
   
    rmse_result <- rmse(y_test_num, pred_var_num)
    mae_result <- mae(y_test_num, pred_var_num)
    
    results[,i] <- c(rmse_result,mae_result)
  }
  
  Meanresultsind <- c(rowMeans(results))
  return(Meanresultsind)
}


result_dfmsvar_lasso <- matrix(NA, 2, length(filenames), dimnames = list(c("rmse","mae"),
                                                                   c(as.list(filenames))))

# iterations for each file
for (i in 1:length(filenames)) {
  filenamesuse <- filenames[i]
  output <- dfmsvar_lasso(filenamesuse) 
  result_dfmsvar_lasso[,i] = output
  
}
mean8 <- rowMeans(result_dfmsvar_lasso)
end.time8 <- Sys.time()
time8 <- end.time8 - start.time8

###########################
# 9. LRPS (lasso)====
start.time9 <- Sys.time()
LRPS_lasso <- function(filename) {
  y <- read.csv(file = filename, header = TRUE)
  y <- data.matrix(y)
  results = array(0, dim = c(2, (nrow(y) - rollingwindow - 1))) # where to store the results
  for(i in 1:(nrow(y) - rollingwindow - 1)) {
    y_train <- y[i: (i-1+rollingwindow),]   #training data
    y_test <- y[(i+ rollingwindow),]
    
    pc <- prcomp(y_train,
                 center = TRUE,
                 scale. = TRUE)
    y_train <- t(t(pc$x[,1:K] %*% t(pc$rotation[,1:K])) * pc$scale + pc$center)
    y_train <- as.data.frame(y_train) #new train 
    # fit sparseVAR
    sparsevar <- fitVAR(y_train, p=1, penalty = "ENET", method = "cv") #default = lasso
    pred_var <- computeForecasts(sparsevar, 1)
    pred_var <- as.data.frame(t(pred_var))
    pred_var <- pred_var[1,]
    
    pred_var_num <- as.numeric(unlist(pred_var))
    
    y_test_num <- as.numeric(unlist(y_test))
    rmse_result <- rmse(y_test_num, pred_var_num)
    mae_result <- mae(y_test_num, pred_var_num)
    results[,i] <- c(rmse_result,mae_result)
  }
  Meanresultsind <- c(rowMeans(results))
  return(Meanresultsind)
}

# Build a container matrix
# Container matrix: rows = MSE,RMSE...; cols = Files 
result_LRPS_lasso <- matrix(NA, 2, length(filenames), dimnames = list(c( "rmse","mae"),
                                                                        c(as.list(filenames))))

# iterations for each file
for (i in 1:length(filenames)) {
  filenamesuse <- filenames[i]
  output <- LRPS_lasso(filenamesuse) 
  result_LRPS_lasso[,i] = output
  
}
mean9 <- rowMeans(result_LRPS_lasso)
end.time9 <- Sys.time()
time9 <- end.time9 - start.time9
time9

################################
# 10. FLsel lasso====
start.time10 <- Sys.time()
FLsel_lasso <- function(filename) {
  y <- read.csv(file = filename, header = TRUE)
  results = array(0, dim = c(2, (nrow(y) - rollingwindow - 1))) # where to store the results
  
  for(i in 1:(nrow(y) - rollingwindow - 1)) {
    y_train <- y[i: (i-1+rollingwindow),]   #training data
    y_test <- y[(i+rollingwindow),]
    
    pc <- prcomp(y_train,
                 center = TRUE,
                 scale. = TRUE)
    y_k <- t(t(pc$x[,1:K] %*% t(pc$rotation[,1:K])) * pc$scale + pc$center)
    y_res <- y_train - y_k
    
    loads <- pc$rotation[,1:K]
    pred_pc <- (loads*y_test)
    
    sparsevar <- fitVAR(y_res, p=1, penalty = "ENET", method = "cv") #default = lasso
    pred_var <- computeForecasts(sparsevar, 1)
    pred_var_value <- t(pred_var)
    pred_var_value <- pred_var + pred_pc
    pred_var_value_num <- as.numeric(unlist(pred_var_value))
    
    y_test_num <- as.numeric(unlist(y_test))
    
    rmse_result <- rmse(y_test_num, pred_var_value_num)
    mae_result <- mae(y_test_num, pred_var_value_num)
    results[,i] <- c(rmse_result,mae_result)
  }
  
  Meanresultsind <- c(rowMeans(results))
  return(Meanresultsind)
}


result_FLsel_lasso<- matrix(NA, 2, length(filenames), dimnames = list(c( "rmse","mae"),
                                                                       c(as.list(filenames))))

# iterations for each file
for (i in 1:length(filenames)) {
  filenamesuse <- filenames[i]
  output <- FLsel_lasso(filenamesuse) 
  result_FLsel_lasso[,i] = output
  
}
mean10 <- rowMeans(result_FLsel_lasso)
end.time10 <- Sys.time()
time10 <- end.time10 - start.time10

#########
mean1
mean2
mean3
mean4
mean5
mean6

mean1
mean7
mean8
mean9
mean10
load("result.RData")


