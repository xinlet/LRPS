library(forecast)
library(dplyr)
library(tidyr)
library(Metrics)
library(matrixStats)
library(VARshrink)
#Data directory=====================================
setwd("~your data directory")
filenames <- list.files(path = "~your data directory",  
                        ("csv"))
K = #num of factors   ##In simulations, K are set. In application studies, K can be obtained by "Choose_k.R"
rollingwindow = #size of the rollingwindow

##################################
LRPS_ridge <- function(filename) {
  y <- read.csv(file = filename, header = TRUE)
  y <- data.matrix(y)
  results = array(0, dim = c(2, (nrow(y) - rollingwindow - 1)))
  for(i in 1:(nrow(y) - rollingwindow - 1)) {
    y_train <- y[i: (i-1+rollingwindow),]  
    y_test <- y[(i+ rollingwindow),]
    
    pc <- prcomp(y_train,
                 center = TRUE,
                 scale. = TRUE)
    y_train <- t(t(pc$x[,1:K] %*% t(pc$rotation[,1:K])) * pc$scale + pc$center)
    y_train <- as.data.frame(y_train) 
    
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

result_LRPS <- matrix(NA, 2, length(filenames), dimnames = list(c( "rmse","mae"),
                                                                        c(as.list(filenames))))

# iterations for each file
for (i in 1:length(filenames)) {
  filenamesuse <- filenames[i]
  output <- LRPS_ridge(filenamesuse) 
  result_LRPS[,i] = output
  
}
mean_over_subject <- rowMeans(result_LRPS)
