library(RSpectra)
library(pls)
library(pryr)

set.seed(123)


# p > q 
lrps_fit <- function(X, Y, k = 1){
  results <- svds(Y, k)
  tmp1 <- solve(crossprod(X))
  tmp2 <- crossprod(X, results$u)
  tmp3 <- tcrossprod(results$d, results$v)
  tmp4 <- (tmp1 %*% tmp2) %*% tmp3
  return(tmp4)
}

RRR_fit <- function(X, Y, k = 1){
  tmp1 <- solve(crossprod(X))
  tmp2 <- crossprod(X, Y)
  Bhat <- tmp1 %*% tmp2
  XBhat <- X %*% Bhat
  svdtmp <- svds(XBhat, k = k)
  tmp3 <- tcrossprod(svdtmp$v)
  Bhat <- Bhat %*% tmp3
  return(Bhat)
}

OLS_fit <- function(X, Y){
  tmp1 <- solve(crossprod(X))
  tmp2 <- crossprod(X, Y)
  Bhat <- tmp1 %*% tmp2
  return(Bhat)
}

# ============================
# NEW METHODS ADDED (PCR, PLS, CCA)
# ============================

PCR_fit <- function(X, Y, k = 1) {
  s <- svds(scale(X, center = TRUE, scale = TRUE), k = k)
  Z <- s$u %*% s$d
  Bz <- solve(crossprod(Z), crossprod(Z, Y))
  Bhat <- s$v %*% Bz
  return(Bhat)
}

PLS_fit <- function(X, Y, k = 1) {
  model <- plsr(Y ~ X, ncomp = k, center = TRUE, scale = TRUE, method = "simpls")
  coef(model, ncomp = k)[,,1]
}

CCA_fit <- function(X, Y, k = 1) {
  if (ncol(Y) >= nrow(Y)) {
    return(matrix(0, ncol(X), ncol(Y)))  # or NA matrix
  }
  cca <- cancor(X, Y)
  A <- cca$xcoef[, 1:k, drop = FALSE]
  B <- cca$ycoef[, 1:k, drop = FALSE]
  Zx <- X %*% A
  Zy <- Y %*% B
  C <- solve(crossprod(Zx), crossprod(Zx, Zy))
  Bhat <- A %*% C %*% t(B)
  return(Bhat)
}

# Main Experiment 

run_experiment <- function(vary = c("n", "p", "q"),
                           base_n = 200, base_p = 20, base_q = 10,
                           scaling_factors = 1:5, reps = 100) {
  
  vary <- match.arg(vary)
  avg_times <- matrix(0, nrow = length(scaling_factors), ncol = 6)
  colnames(avg_times) <- c("OLS", "RRR", "LRPS", "PCR", "PLS", "CCA")
  rownames(avg_times) <- scaling_factors
  
  all_raw_times <- list()
  all_raw_mem   <- list()
  
  for (i_idx in seq_along(scaling_factors)) {
    i <- scaling_factors[i_idx]
    n <- base_n; p <- base_p; q <- base_q
    if (vary == "n") n <- base_n * i
    if (vary == "p") p <- base_p * i
    if (vary == "q") q <- base_q * i
    
    times <- matrix(0, nrow = reps, ncol = 6)
    colnames(times) <- c("OLS", "RRR", "LRPS", "PCR", "PLS", "CCA")
    
    mem_use <- matrix(0, nrow = reps, ncol = 6)
    colnames(mem_use) <- c("OLS", "RRR", "LRPS", "PCR", "PLS", "CCA")
    
    for (r in 1:reps) {
      X <- matrix(rnorm(n * p), n, p)
      Y <- matrix(rnorm(n * q), n, q)
      
      # ---- OLS ----
      mem_use[r, "OLS"]  <- mem_change(tmp <- OLS_fit(X, Y))
      times[r, "OLS"]    <- system.time(OLS_fit(X, Y))[3]
      
      # ---- RRR ----
      mem_use[r, "RRR"]  <- mem_change(tmp <- RRR_fit(X, Y))
      times[r, "RRR"]    <- system.time(RRR_fit(X, Y))[3]
      
      # ---- LRPS ----
      mem_use[r, "LRPS"] <- mem_change(tmp <- lrps_fit(X, Y))
      times[r, "LRPS"]   <- system.time(lrps_fit(X, Y))[3]
      
      # ---- PCR ----
      mem_use[r, "PCR"]  <- mem_change(tmp <- PCR_fit(X, Y, k = 1))
      times[r, "PCR"]    <- system.time(PCR_fit(X, Y, k = 1))[3]
      
      # ---- PLS ----
      mem_use[r, "PLS"]  <- mem_change(tmp <- PLS_fit(X, Y, k = 1))
      times[r, "PLS"]    <- system.time(PLS_fit(X, Y, k = 1))[3]
      
      # ---- CCA ----
      mem_use[r, "CCA"]  <- mem_change(tmp <- CCA_fit(X, Y, k = 1))
      times[r, "CCA"]    <- system.time(CCA_fit(X, Y, k = 1))[3]
    }
    
    avg_times[i_idx, ] <- apply(times, 2, mean)
    all_raw_times[[i_idx]] <- times
    all_raw_mem[[i_idx]]   <- mem_use
  }
  
  return(list(avg_times = avg_times,
              raw_times = all_raw_times,
              raw_memory = all_raw_mem))
}


# Run Experiments (UNCHANGED)

result_n <- run_experiment("n")
result_p <- run_experiment("p")

get_relative_times <- function(avg_times) {
  avg_times / avg_times[, "OLS"]
}

rel_n <- get_relative_times(result_n$avg_times)
rel_p <- get_relative_times(result_p$avg_times)


# Plotting 

plot_relative_time <- function(rel_times, title_text, y_limits = NULL) {
  scaling_factors <- 1:nrow(rel_times)
  
  if (is.null(y_limits)) {
    y_limits <- range(rel_times)
  }
  
  plot(scaling_factors, rel_times[, 1], type = "b", col = "blue", pch = 16,
       ylim = y_limits, xlab = "Scaling Factor",
       ylab = "Relative Time", main = title_text)
  
  for (j in 2:ncol(rel_times)) {
    lines(scaling_factors, rel_times[, j], type = "b",
          pch = 15 + j, col = j)
  }
  
  legend("topright", legend = colnames(rel_times),
         col = 1:ncol(rel_times), pch = 16:(15 + ncol(rel_times)), lty = 1)
}

par(mfrow = c(1, 2))
plot_relative_time(rel_n, "Increasing n", y_limits = c(0, 7))
plot_relative_time(rel_p, "Increasing p", y_limits = c(0, 7))


pdf(file = "complex_fig1.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
plot_relative_time(rel_n, "Increasing n", y_limits = c(0, 7))
plot_relative_time(rel_p, "Increasing p", y_limits = c(0, 7))
dev.off()




# Self-relative times 

get_self_relative_times <- function(avg_times) {
  sweep(avg_times, 2, avg_times[1, ], FUN = "/")
}

rel_n_self <- get_self_relative_times(result_n$avg_times)
rel_p_self <- get_self_relative_times(result_p$avg_times)

plot_self_relative_time <- function(rel_times, title_text, y_limits = NULL) {
  scaling_factors <- 1:nrow(rel_times)
  
  if (is.null(y_limits)) {
    y_limits <- range(rel_times, finite = TRUE)
  }
  
  # Colors and symbols in correct order
  cols <- c("blue", 2:(ncol(rel_times)))   # first is blue, rest are 2,3,...
  pchs <- c(16, 17:(15 + ncol(rel_times) - 1))
  
  # First line
  plot(scaling_factors, rel_times[, 1], type = "b",
       col = cols[1], pch = pchs[1],
       ylim = y_limits, xlab = "Scaling Factor",
       ylab = "Relative Time (vs. unit scale)", main = title_text)
  
  # Remaining lines
  for (j in 2:ncol(rel_times)) {
    lines(scaling_factors, rel_times[, j], type = "b",
          col = cols[j], pch = pchs[j])
  }
  
  legend("topleft", legend = colnames(rel_times),
         col = cols, pch = pchs, lty = 1)
}


par(mfrow = c(1, 2))
plot_self_relative_time(rel_n_self, "Increasing n ", y_limits = c(0, 12))
plot_self_relative_time(rel_p_self, "Increasing p ", y_limits = c(0, 12))

pdf(file = "complex_fig2.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
plot_self_relative_time(rel_n_self, "Increasing n ", y_limits = c(0, 12))
plot_self_relative_time(rel_p_self, "Increasing p ", y_limits = c(0, 12))
dev.off()


# Memory Usage Plot 

# Get average memory across reps for each scale
get_avg_memory <- function(raw_memory_list, to_MB = TRUE) {
  k <- length(raw_memory_list)
  m <- ncol(raw_memory_list[[1]])
  
  avg_mem <- matrix(0, nrow = k, ncol = m)
  colnames(avg_mem) <- colnames(raw_memory_list[[1]])
  rownames(avg_mem) <- paste0("scale_", 1:k)
  
  for (i in 1:k) {
    avg_mem[i, ] <- colMeans(raw_memory_list[[i]], na.rm = TRUE)
  }
  
  if (to_MB) {
    avg_mem <- avg_mem / (1024^2)   # convert bytes → MB
  }
  
  return(avg_mem)
}
get_avg_times <- function(avg_times_matrix) {
  out <- avg_times_matrix
  rownames(out) <- paste0("scale_", 1:nrow(avg_times_matrix))
  return(out)
}

make_time_memory_table <- function(result_obj) {
  
  avg_time <- get_avg_times(result_obj$avg_times)
  avg_mem  <- get_avg_memory(result_obj$raw_memory)
  
  methods <- colnames(avg_time)
  scales  <- rownames(avg_time)
  
  df_list <- list()
  
  for (s in seq_along(scales)) {
    for (m in seq_along(methods)) {
      df_list[[length(df_list) + 1]] <- data.frame(
        Scale = scales[s],
        Method = methods[m],
        Avg_Time_sec = avg_time[s, m],
        Avg_Memory_MB = avg_mem[s, m]
      )
    }
  }
  
  return(do.call(rbind, df_list))
}

# Fix memory for a result object
fix_memory <- function(result_obj) {
  for (i in seq_along(result_obj$raw_memory)) {
    result_obj$raw_memory[[i]] <- abs(result_obj$raw_memory[[i]])
  }
  return(result_obj)
}

result_n_fixed <- fix_memory(result_n)
result_p_fixed <- fix_memory(result_p)

table_n <- make_time_memory_table(result_n_fixed)
table_p <- make_time_memory_table(result_p_fixed)

table_n
table_p




#q > p====
set.seed(123)
lrps_fit <- function(X, Y, k = 1){
  results <- svds(Y, k)
  tmp1 <- solve(crossprod(X))
  tmp2 <- crossprod(X, results$u)
  tmp3 <- tmp1 %*% tmp2
  tmp4 <- tmp3 %*% results$d
  tmp5 <- tcrossprod(tmp4, results$v)
  return(tmp5)
}


run_experiment <- function(vary = c("n", "p", "q"),
                           base_n = 60, base_p = 30, base_q = 100,
                           scaling_factors = 1:5, reps = 100) {
  
  vary <- match.arg(vary)
  avg_times <- matrix(0, nrow = length(scaling_factors), ncol = 6)
  colnames(avg_times) <- c("OLS", "RRR", "LRPS", "PCR", "PLS", "CCA")
  rownames(avg_times) <- scaling_factors
  
  all_raw_times <- list()
  all_raw_mem   <- list()
  
  for (i_idx in seq_along(scaling_factors)) {
    i <- scaling_factors[i_idx]
    n <- base_n; p <- base_p; q <- base_q
    if (vary == "n") n <- base_n * i
    if (vary == "p") p <- base_p * i
    if (vary == "q") q <- base_q * i
    
    times <- matrix(0, nrow = reps, ncol = 6)
    colnames(times) <- c("OLS", "RRR", "LRPS", "PCR", "PLS", "CCA")
    
    mem_use <- matrix(0, nrow = reps, ncol = 6)
    colnames(mem_use) <- c("OLS", "RRR", "LRPS", "PCR", "PLS", "CCA")
    
    for (r in 1:reps) {
      X <- matrix(rnorm(n * p), n, p)
      Y <- matrix(rnorm(n * q), n, q)
      
      # ---- OLS ----
      mem_use[r, "OLS"]  <- mem_change(tmp <- OLS_fit(X, Y))
      times[r, "OLS"]    <- system.time(OLS_fit(X, Y))[3]
      
      # ---- RRR ----
      mem_use[r, "RRR"]  <- mem_change(tmp <- RRR_fit(X, Y))
      times[r, "RRR"]    <- system.time(RRR_fit(X, Y))[3]
      
      # ---- LRPS ----
      mem_use[r, "LRPS"] <- mem_change(tmp <- lrps_fit(X, Y))
      times[r, "LRPS"]   <- system.time(lrps_fit(X, Y))[3]
      
      # ---- PCR ----
      mem_use[r, "PCR"]  <- mem_change(tmp <- PCR_fit(X, Y, k = 1))
      times[r, "PCR"]    <- system.time(PCR_fit(X, Y, k = 1))[3]
      
      # ---- PLS ----
      mem_use[r, "PLS"]  <- mem_change(tmp <- PLS_fit(X, Y, k = 1))
      times[r, "PLS"]    <- system.time(PLS_fit(X, Y, k = 1))[3]
      
      # ---- CCA ----
      mem_use[r, "CCA"]  <- mem_change(tmp <- CCA_fit(X, Y, k = 1))
      times[r, "CCA"]    <- system.time(CCA_fit(X, Y, k = 1))[3]
    }
    
    avg_times[i_idx, ] <- apply(times, 2, mean)
    all_raw_times[[i_idx]] <- times
    all_raw_mem[[i_idx]]   <- mem_use
  }
  
  return(list(avg_times = avg_times,
              raw_times = all_raw_times,
              raw_memory = all_raw_mem))
}


get_relative_times <- function(avg_times) {
  avg_times / avg_times[, "OLS"]
}
result_n <- run_experiment("n")
result_q <- run_experiment("q")


rel_n <- get_relative_times(result_n$avg_times)
rel_q <- get_relative_times(result_q$avg_times)
plot_relative_time <- function(rel_times, title_text, y_limits = NULL) {
  scaling_factors <- 1:nrow(rel_times)
  
  if (is.null(y_limits)) {
    y_limits <- range(rel_times)
  }
  
  plot(scaling_factors, rel_times[, 1], type = "b", col = "blue", pch = 16,
       ylim = y_limits, xlab = "Scaling Factor",
       ylab = "Relative Time", main = title_text)
  
  for (j in 2:ncol(rel_times)) {
    lines(scaling_factors, rel_times[, j], type = "b",
          pch = 15 + j, col = j)
  }
  
  legend("topleft", legend = colnames(rel_times),
         col = 1:ncol(rel_times), pch = 16:(15 + ncol(rel_times)), lty = 1)
}

par(mfrow = c(1, 2))
plot_relative_time(rel_n, "Increasing n", y_limits = c(0, 20))
plot_relative_time(rel_q, "Increasing q", y_limits = c(0, 20))

par(mfrow = c(1, 2))
plot_relative_time(rel_n[, colnames(rel_n) != "CCA"], 
                        "Increasing n", 
                        y_limits = c(0, 10))
plot_relative_time(rel_q[, colnames(rel_q) != "CCA"], 
                        "Increasing q", 
                        y_limits = c(0, 10))

pdf(file = "complex_fig3.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
plot_relative_time(rel_n[, colnames(rel_n) != "CCA"], 
                   "Increasing n", 
                   y_limits = c(0, 12))
plot_relative_time(rel_q[, colnames(rel_q) != "CCA"], 
                   "Increasing q", 
                   y_limits = c(0, 12))
dev.off()


# Self-relative times 

get_self_relative_times <- function(avg_times) {
  sweep(avg_times, 2, avg_times[1, ], FUN = "/")
}

rel_n_self <- get_self_relative_times(result_n$avg_times)
rel_q_self <- get_self_relative_times(result_q$avg_times)

plot_self_relative_time <- function(rel_times, title_text, y_limits = NULL) {
  scaling_factors <- 1:nrow(rel_times)
  
  if (is.null(y_limits)) {
    y_limits <- range(rel_times, finite = TRUE)
  }
  
  # Colors and symbols in correct order
  cols <- c("blue", 2:(ncol(rel_times)))   # first is blue, rest are 2,3,...
  pchs <- c(16, 17:(15 + ncol(rel_times) - 1))
  
  # First line
  plot(scaling_factors, rel_times[, 1], type = "b",
       col = cols[1], pch = pchs[1],
       ylim = y_limits, xlab = "Scaling Factor",
       ylab = "Relative Time (vs. unit scale)", main = title_text)
  
  # Remaining lines
  for (j in 2:ncol(rel_times)) {
    lines(scaling_factors, rel_times[, j], type = "b",
          col = cols[j], pch = pchs[j])
  }
  
  legend("topright", legend = colnames(rel_times),
         col = cols, pch = pchs, lty = 1)
}



par(mfrow = c(1, 2))
plot_self_relative_time(rel_n_self, "Increasing n ", y_limits = c(0, 12))
plot_self_relative_time(rel_q_self, "Increasing q ", y_limits = c(0, 12))

par(mfrow = c(1, 2))
plot_self_relative_time(rel_n_self[, colnames(rel_n_self) != "CCA"], 
                        "Increasing n", 
                        y_limits = c(0, 10))
plot_self_relative_time(rel_q_self[, colnames(rel_q_self) != "CCA"], 
                        "Increasing q", 
                        y_limits = c(0, 10))



pdf(file = "complex_fig4.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
plot_self_relative_time(rel_n_self[, colnames(rel_n_self) != "CCA"], 
                        "Increasing n", 
                        y_limits = c(0, 10))
plot_self_relative_time(rel_q_self[, colnames(rel_q_self) != "CCA"], 
                        "Increasing q", 
                        y_limits = c(0, 10))
dev.off()


# Memory Usage Plot
# Get average memory across reps for each scale
get_avg_memory <- function(raw_memory_list, to_MB = TRUE) {
  k <- length(raw_memory_list)
  m <- ncol(raw_memory_list[[1]])
  
  avg_mem <- matrix(0, nrow = k, ncol = m)
  colnames(avg_mem) <- colnames(raw_memory_list[[1]])
  rownames(avg_mem) <- paste0("scale_", 1:k)
  
  for (i in 1:k) {
    avg_mem[i, ] <- colMeans(raw_memory_list[[i]], na.rm = TRUE)
  }
  
  if (to_MB) {
    avg_mem <- avg_mem / (1024^2)   # convert bytes → MB
  }
  
  return(avg_mem)
}
get_avg_times <- function(avg_times_matrix) {
  out <- avg_times_matrix
  rownames(out) <- paste0("scale_", 1:nrow(avg_times_matrix))
  return(out)
}

make_time_memory_table <- function(result_obj) {
  
  avg_time <- get_avg_times(result_obj$avg_times)
  avg_mem  <- get_avg_memory(result_obj$raw_memory)
  
  methods <- colnames(avg_time)
  scales  <- rownames(avg_time)
  
  df_list <- list()
  
  for (s in seq_along(scales)) {
    for (m in seq_along(methods)) {
      df_list[[length(df_list) + 1]] <- data.frame(
        Scale = scales[s],
        Method = methods[m],
        Avg_Time_sec = avg_time[s, m],
        Avg_Memory_MB = avg_mem[s, m]
      )
    }
  }
  
  return(do.call(rbind, df_list))
}

# Fix memory for a result object
fix_memory <- function(result_obj) {
  for (i in seq_along(result_obj$raw_memory)) {
    result_obj$raw_memory[[i]] <- abs(result_obj$raw_memory[[i]])
  }
  return(result_obj)
}

result_n_fixed <- fix_memory(result_n)
result_q_fixed <- fix_memory(result_q)

table_n <- make_time_memory_table(result_n_fixed)
table_q <- make_time_memory_table(result_q_fixed)

table_n
table_q

