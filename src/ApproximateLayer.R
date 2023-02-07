# Flush r-studio and set working directory 
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
set.seed(0)

# read in files

MSA = read.csv("UsedMSANS1.csv", header = F, colClasses = "character")
layer_values = read.csv("NodeValues2L3N.csv", header = F)
num_layers = dim(layer_values)[1]
num_seq = dim(layer_values)[2]
assertthat::are_equal(num_seq, dim(MSA)[1])

# create approximation for first layer (amino acid to node values)
# create one data frame per num_layers
library(R.ROSETTA)
library(arules)
library(varhandle)

lay1 <- list()
for (i in 1:num_layers) {
  dec = discretize(unlist(unname(layer_values[i,])), method = "frequency", breaks = 3, labels = c(1,2,3))
  df <- data.frame(cbind(MSA[,2:250], unfactor(dec)))
  print(df[1:3,])
  lay1[[i]] <- rosetta(df, discrete = TRUE, clroc = c(1,2))
}

head(viewRules(lay1[[1]]$main[lay1[[1]]$main$pValue<= 0.05,]))


val <- apply(unname(layer_values), 1, discretize, method = "frequency", breaks = 3, labels = c(1,2,3))
dt <- data.frame(val, as.numeric(MSA[,251]))
lay2 <- rosetta(dt, discrete = TRUE, clroc = c(-1, 1))
