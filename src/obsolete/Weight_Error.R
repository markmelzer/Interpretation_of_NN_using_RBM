# Flush r-studio and set working directory
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
set.seed(0)


# load necessary packages
library(R.ROSETTA)
library(arules)


# load data
data = read.csv("../Neuron1_changes.csv", header = F)


# discretize
discretize_data <- function(data, params = list(method = "interval", breaks = 2, 
                                                labels = c(0, 1))){
  return(discretizeDF(data, default = params))
}

df = discretize_data(data)

df = data
df[df > 0] = 1 
df[df < 0] = -1 


# run rosetta

ros <- rosetta(df, discrete = T, roc = T, clroc = 1)
