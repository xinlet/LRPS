library(tidyr)
library(ggplot2)
library(ggeasy)
library(ggpubr)

#Data directory=====================================
setwd("~Atlanta data directory") #change to your own directory
filenames <- list.files(path = "~Atlanta data directory", ("csv")) # change to your own directory
filenames <- filenames[1:length(filenames)]
rollingwindow = 120 #in our study is 120
#====
data <- read.csv(filenames[1], header = TRUE) 
mindim <- min(rollingwindow, ncol(data))

analyzek <- function(filename) {
  y <- read.csv(file = filename, header = TRUE)
  y <- data.matrix(y)
  results = array(0, dim = c(mindim, (nrow(y) - rollingwindow - 1))) 
  for(i in 1:(nrow(y) - rollingwindow - 1)) {
    y_train <- y[i: (i-1+rollingwindow),]   
    
    pc <- prcomp(y_train,
                 center = TRUE,
                 scale. = TRUE)
    choose_k <- cumsum(pc$sdev^2 / sum(pc$sdev^2))
    results[,i] <- c(choose_k)
  }
  Meanresultsind <- c(rowMeans(results))
  return(Meanresultsind)
}

result_k  <- matrix(NA, mindim, length(filenames))

for (i in 1:length(filenames)) {
  filenamesuse <- filenames[i]
  output <- analyzek(filenamesuse) 
  result_k[,i] = output
}

k_mean <- rowMeans(result_k)
k_num <- c(1:mindim)
matrix_k <- as.data.frame(cbind(k_num, k_mean))

ggplot(matrix_k) + geom_line(aes(x = k_num, y = k_mean)) +
  labs(y= "Cumulative explained Variance", x = "Number of k") +
  scale_x_continuous(breaks = seq(0, mindim, by = 10)) +
  geom_hline(yintercept=0.9, linetype="dashed", color = "red") + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) 

K <- which(matrix_k[,2] > 0.9)[1]

