library(tidyr)
library(ggplot2)
library(ggeasy)
library(ggpubr)
#Data directory=====================================
setwd("~/your data directory")
filenames <- list.files(path = "~/your data directory",  
                        ("csv"))
filenames <- filenames[1:length(filenames)]
rollingwindow = 50 #your defined rolling window size, in our study is either 50 or 120
#====

analyzek <- function(filename) {
  y <- read.csv(file = filename, header = TRUE)
  y <- data.matrix(y)
  results = array(0, dim = c( min(ncol(y), rollingwindow), (nrow(y) - rollingwindow - 1)))
  for(i in 1:(nrow(y) - rollingwindow - 1)) {
    y_train <- y[i: (i-1+rollingwindow),]   #training data
    #y_test <- y[(i+ rollingwindow),]
    
    pc <- prcomp(y_train,
                 center = TRUE,
                 scale. = TRUE)
    choose_k <- cumsum(pc$sdev^2 / sum(pc$sdev^2))
    results[,i] <- c(choose_k)
  }
  Meanresultsind <- c(rowMeans(results))
  return(Meanresultsind)
}

result_k  <- matrix(NA, min(ncol(y), rollingwindow), length(filenames))

for (i in 1:length(filenames)) {
  filenamesuse <- filenames[i]
  output <- analyzek(filenamesuse) 
  result_k[,i] = output
}

k_mean <- rowMeans(result_k)
k_num <- c(1:min(ncol(y), rollingwindow))
matrix_k <- as.data.frame(cbind(k_num, k_mean))

ggplot(matrix_k) + geom_line(aes(x = k_num, y = k_mean)) +
  labs(y= "Cumulative explained Variance", x = "Number of k") +
  scale_x_continuous(breaks = seq(0, min(ncol(y), rollingwindow), by = 10)) +
  geom_hline(yintercept=0.9, linetype="dashed", color = "red") + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) 

