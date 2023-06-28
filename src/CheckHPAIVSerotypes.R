# Flush r-studio and set working directory 
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
set.seed(0)

# load the necessary packages
library(dplyr)

# read the MSA
df <- read.csv("../data/NS1/NS1.csv", header = T, colClasses = "character")

# extract serotype of objects
serotype <- unlist(lapply(df$id, function(x){unlist(strsplit(x, '[*]'))[5]}))

# extract HA serotype
H_serotype <- unlist(lapply(serotype, function(x){unlist(strsplit(x, 'N'))[1]}))

H_serotype %>% table

# extract objects of serotype H5 or H7
idx <- which(H_serotype %in% c("H5", "H7"))
dat <- df[idx,]
dim(dat)

# save data
write.csv(dat, "../data/NS1/NS1_H5_H7.csv", quote = F, row.names = F)
