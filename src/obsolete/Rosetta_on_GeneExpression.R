# Flush r-studio and set working directory 
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
set.seed(0)


# get data
data = read.csv("../GE_train.csv", header = F)



# standardize data
features <- data[,1:100]
boxplot(features)


standardize <- function(x){
  return((x - mean(x))/(sd(x)))
}

features <- apply(features, 2, standardize)
boxplot(features)

#boxplot(t(features))
#features <- apply(features, 1, standardize)
#boxplot(t(features))

# discretization
library(arules)
feat_discrete <- data.frame(apply(data[,1:100], 2, discretize, method = "frequency",
                       breaks = 3, labels = c("1", "2", "3")))



df_discrete <- cbind(feat_discrete, data[,101])
names(df_discrete) <- append(names(feat_discrete), "decision")




# Rosetta
library(R.ROSETTA)
out1 <- rosetta(df_discrete, discrete = TRUE, reducer = "Johnson")
