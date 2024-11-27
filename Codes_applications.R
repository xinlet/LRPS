#load("~.rda")

#For UK and Beijing datasets
split_index <- floor(0.8 * nrow(X))
Xtrain <- X[1:split_index, ]
Ytrain <- Y[1:split_index, ]
Xtest <- X[(split_index + 1):nrow(X), ]
Ytest <- Y[(split_index + 1):nrow(Y), ]

#For USA dataset
Xtrain <- X[1:578, ]
Ytrain <- Y[1:578, ]
Xtest <- X[(578 + 1):nrow(X), ]
Ytest <- Y[(578 + 1):nrow(Y), ]

#For Gene dataset
Xtrain <- X[1:100, ]
Ytrain <- Y[1:100, ]
Xtest <- X[(100 + 1):nrow(X), ]
Ytest <- Y[(100 + 1):nrow(Y), ]

#OLS computing====
# Compute OLS estimates on training data
tmp <- solve(crossprod(Xtrain)) %*% t(Xtrain)
Bhat_ols <- tmp %*% Ytrain
Yhat <- Xtest %*% Bhat_ols
nt <- nrow(Ytest)
q <- ncol(Ytest)
prediction_error_OLS <- sum((Yhat - Ytest)^2) / (nt * q)
cat("Prediction error for OLS:", prediction_error_OLS)

#RRR computing====
getMLEs <- function(X, Y) {
  if (!is.matrix(X)) {
    X <- matrix(X)
  }
  n <- nrow(Y)
  m <- ncol(X)
  tmp <- solve(crossprod(X)) %*% t(X)
  Bhat <- tmp %*% Y
  Shat <- t(Y) %*% (diag(n) - X %*% tmp) %*% Y / n
  return(list(Bhat = Bhat, Shat = Shat))
}
cross_validate_k_rrr_two_fold <- function(Ytrain, Xtrain) {
  n <- nrow(Ytrain)
  p <- ncol(Xtrain)
  k_values <- 1:min(n, p)  
  set.seed(123)
  half_index <- floor(n / 2)
  indices <- sample(1:n)
  first_half <- indices[1:half_index]
  second_half <- indices[(half_index + 1):n]
  
  mse_k <- numeric(length(k_values))

  for (fold in 1:2) {
    if (fold == 1) {
      train_indices <- first_half
      val_indices <- second_half
    } else {
      train_indices <- second_half
      val_indices <- first_half
    }
    
    Ytrain_fold <- Ytrain[train_indices, , drop = FALSE]
    Xtrain_fold <- Xtrain[train_indices, , drop = FALSE]
    Yval_fold <- Ytrain[val_indices, , drop = FALSE]
    Xval_fold <- Xtrain[val_indices, , drop = FALSE]
    
    for (k_index in 1:length(k_values)) {
      k <- k_values[k_index]
      
      MLEfull <- getMLEs(X = Xtrain_fold, Y = Ytrain_fold)
      
      tmp <- Xtrain_fold %*% MLEfull$Bhat
      svdtmp <- svd(tmp)
      vk <- svdtmp$v[, 1:k]
      Bhat_reduced <- MLEfull$Bhat %*% vk %*% t(vk)
      
      Yhat_val <- Xval_fold %*% Bhat_reduced
      
      mse_k[k_index] <- mse_k[k_index] + mean((Yval_fold - Yhat_val)^2) / 2  # Average over 2 folds
    }
  }
  
  optimal_k <- k_values[which.min(mse_k)]
  
  return(list(optimal_k = optimal_k, avg_mse_k = mse_k))
}

optimal_k_rrr = cross_validate_k_rrr_two_fold(Ytrain, Xtrain)$optimal_k
MLEfull = getMLEs(X = Xtrain, Y = Ytrain)
tmp = Xtrain %*% MLEfull$Bhat
svdtmp = svd(tmp)
tmp2 = svdtmp$v
vk = tmp2[,1:optimal_k_rrr]
Bhatf = MLEfull$Bhat %*% vk %*% t(vk)
Yhat <- Xtest %*% Bhatf

# Compute prediction error
nt <- nrow(Ytest)
q <- ncol(Ytest)
prediction_error_rrr <- sum((Yhat - Ytest)^2) / (nt * q)
prediction_error_rrr
cat("Prediction error for RRR:", prediction_error_rrr)

#LRPS computing====
cross_validate_k_two_fold <- function(Ytrain, Xtrain) {
  n <- nrow(Ytrain)
  q <- ncol(Ytrain)
  k_values <- 1:min(n/2, q)  
  set.seed(123)
  half_index <- floor(n / 2)
  indices <- sample(1:n)  
  first_half <- indices[1:half_index]
  second_half <- indices[(half_index + 1):n]
  
  mse_matrix <- matrix(NA, nrow = length(k_values), ncol = 2)
  
  for (fold in 1:2) {
    if (fold == 1) {
      train_indices <- first_half
      test_indices <- second_half
    } else {
      train_indices <- second_half
      test_indices <- first_half
    }
    
    Ytrain_fold <- Ytrain[train_indices, , drop = FALSE]
    Xtrain_fold <- Xtrain[train_indices, , drop = FALSE]
    Ytest_fold <- Ytrain[test_indices, , drop = FALSE]
    Xtest_fold <- Xtrain[test_indices, , drop = FALSE]
    
    mse_k <- numeric(length(k_values))
    
    for (k_index in 1:length(k_values)) {
      k <- k_values[k_index]
      pc <- prcomp(Ytrain_fold, center = FALSE, scale. = FALSE)
      Ys_fold <- t(t(pc$x[, 1:k] %*% t(pc$rotation[, 1:k])))
      tmp <- solve(crossprod(Xtrain_fold)) %*% t(Xtrain_fold)
      Bhat <- tmp %*% Ys_fold
      Yhat_test <- Xtest_fold %*% Bhat
      mse_k[k_index] <- mean(as.vector((Ytest_fold - Yhat_test)^2))
    }
    
    mse_matrix[, fold] <- mse_k
  }
  
  avg_mse_k <- rowMeans(mse_matrix)
  
  optimal_k <- k_values[which.min(avg_mse_k)]
  
  return(list(optimal_k = optimal_k, avg_mse_k = avg_mse_k, mse_matrix = mse_matrix))
}

optimal_k_lrps = cross_validate_k_two_fold(Ytrain,Xtrain)$optimal_k

pc <- prcomp(Ytrain,
             center = FALSE,
             scale. = FALSE)
Ys <- t(t(pc$x[,1:optimal_k_lrps] %*% t(pc$rotation[,1:optimal_k_lrps])))
tmp <- solve(crossprod(Xtrain)) %*% t(Xtrain)
Bhat <- tmp %*% Ys
Yhat <- Xtest %*% Bhat

# Compute prediction error
nt <- nrow(Ytest)
q <- ncol(Ytest)
prediction_error_lrps <- sum((Yhat - Ytest)^2) / (nt * q)
cat("Prediction error for lrps:", prediction_error_lrps)
