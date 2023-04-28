# Flush r-studio and set working directory 
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
set.seed(0)

# load necessary packages
library(R.ROSETTA)

# The result is saved in ../results/Exp1_HPAIV/InputSizeAnalysis.rds

# get data, remove id's
dat <- read.csv("../data/NS1/NS1.csv", header = T, colClasses = "character")
dat <- dat[, 2:dim(dat)[2]]

# separate positive and negative objects 
pos <- dat[dat$Pathogenicity == "1",]
neg <- dat[dat$Pathogenicity == "0",]

# decide on data split sizes (e.g., 20 splits)
splits <- 20

# arrays to save the result
pos_feat <- c()
neg_feat <- c()
time <- c()
num_obj <- c()

# loop to run Rosetta with different sized inputs
for (i in 1:splits) {
  # get a data set of pre-defined size
  df <- rbind(pos[1:round(i * dim(pos)[1]/splits),], neg[1:round(i * dim(neg)[1]/splits),])
  df <- df[sample(1:dim(df)[1]),]
  num_obj <- append(num_obj, dim(df)[1])
  
  # run Rosetta, recalculate statistics, and extract significant rules
  print(paste0("Run Rosetta for round: ", as.character(i)))
  t0 <- Sys.time()
  ros <- rosetta(df, discrete = T, underSample = T)
  t1 <- Sys.time()
  time <- append(time, as.numeric(t1-t0))
  rec <- recalculateRules(df, ros$main, discrete = T)
  sig <- rec[rec$pValue <= 0.05,]
  
  # get Features for significant rules
  feat <- getFeatures(sig)
  pos_found <- feat$features$"1"
  neg_found <- feat$features$"0"
  
  # append length to result array
  pos_feat <- append(pos_feat, length(pos_found))
  neg_feat <- append(neg_feat, length(neg_found))
  
  # shuffle the data to randomize used objects
  pos <- pos[sample(1:dim(pos)[1]),]
  neg <- neg[sample(1:dim(neg)[1]),]
}

# scale the time data points measured in seconds to minutes
time[1:4] <- time[1:4]/60

# plot the results
library(ggplot2)
library(ggpubr)

# plot the computation time over the number of input objects
ggplot(data = data.frame(num_obj, time), aes(num_obj, time)) +
  geom_point() +
  ylab("Computation time [min]") +
  xlab("Input data size") +
  ggtitle("Computation time over input data size")

# plot the logarithmic computation time over the number of input objects
# add a fitted linear model
ggplot(data = data.frame(num_obj, log(time)), aes(num_obj, log(time))) +
  geom_point() +
  geom_smooth(method = lm) +
  stat_cor(method = "pearson", label.x = 200) +
  ylab("Logarithmic computation time [min]")  +
  xlab("Input data size") +
  ggtitle("Logarithmic computation time over input data size")

# plot the number of found features for HP over the number of input objects
# add a fitted linear model
ggplot(data = data.frame(num_obj, pos_feat), aes(num_obj, pos_feat)) +
  geom_point() +
  geom_smooth(method = lm) +
  stat_cor(method = "pearson", label.x = 200) +
  ylab("Number of features")  +
  xlab("Input data size") +
  ggtitle("Number of features appearing in significant rules for HP over input data size")

# plot the number of found features for LP over the number of input objects
# add a fitted linear model
ggplot(data = data.frame(num_obj, neg_feat), aes(num_obj, neg_feat)) +
  geom_point() +
  geom_smooth(method = lm) +
  stat_cor(method = "pearson", label.x = 200) +
  ylab("Number of features") +
  xlab("Input data size") +
  ggtitle("Number of features appearing in significant rules for LP over input data size")

# save results to RDS
tmp <- data.frame(num_obj, time, pos_feat, neg_feat)
saveRDS(tmp, "../results/Exp1_HPAIV/InputSizeAnalysis.rds")

