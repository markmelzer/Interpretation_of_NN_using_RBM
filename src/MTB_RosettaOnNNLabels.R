# Flush r-studio and set working directory 
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
set.seed(0)

# load necessary packages
library(dplyr)
library(R.ROSETTA)

# ======================= Build RBM on NN output ==============================

# load data
dat <- read.csv("../data/Tubercolosis/NewTrain2.csv", header = F, colClasses = "character")

# change names of data frame to positions and decisions
feat <- readRDS("../data/Tubercolosis/Relevant_attributes.RDS")

names(dat)[1:(dim(dat)[2] - 2)] = feat
names(dat)[dim(dat)[2] - 1] = "True"
names(dat)[dim(dat)[2]] = "Prediction"

# extract true labels and remove from data frame
truth <- dat[,(dim(dat)[2] - 1)]
df <- dat[,c(1:(dim(dat)[2] - 2), dim(dat)[2])]

# run Rosetta on NN output
ros_g <- rosetta(df, discrete = T, underSample = T, reducer = "Genetic")

# get accurate rule parameters
rec <- recalculateRules(df, ros_g$main, discrete = T)
rec <- rec[rec$pValue<=0.05,]

# check rules
viewRules(rec)


# ======================== Test RBM on NN output ==============================

# test Rosetta classifier (rec) on test set
test <- read.csv("../data/Tubercolosis/NewTest.csv", colClasses = "character", header = F)

# rename features to their original name
test_truth <- test[,(dim(test)[2] - 1)]
test_df <- test[,c(1:(dim(test)[2] - 2), dim(test)[2])]
names(test_df) <- c(feat, "Prediction")

pred <- predictClass(test_df[,1:(dim(test_df)[2] - 1)], rec, discrete = T, validate = T, defClass = test_df[,dim(test_df)[2]], normalizeMethod = "rulnum")
pred$accuracy
table(pred$out[,c("currentClass", "predictedClass")])


# ================= Analysis of wrongly classified objects ====================

# extract wrongly classified objects
wrongObj <- which(df[,dim(df)[2]] != truth)

# get supportSet as list
supportSet <- lapply(rec$supportSetRHS, function(x){as.numeric(unlist(strsplit(x, ",")))})

# create heatmap to identify misclassified objects in support set of rules
rules <- list()

for(i in wrongObj){
  tmp = which(unlist(lapply(supportSet, function(x){(i %in% x)})))
  rules[[length(rules) + 1]] <- as.numeric(row.names(rec[tmp,]))
}

r <- unlist(rules) %>% unique

x = lapply(rules, function(x){as.numeric(r %in% x)}) %>% data.frame
names(x) = wrongObj
row.names(x) <- r


heatmap(as.matrix(x), scale = "none", Colv = F, col = c("white", "black"),
        xlab = "Object Number", ylab = "Rule Rank")


# based on heatmap inspect rules 
viewRules(rec[c(1, 2),])

# extract all objects with those objects
df_tmp1 <- filter(dat, V594 == "S")
df_tmp2 <- filter(dat, V594 == "T")

df_tmp1[,c("True", "Prediction")] %>% table
df_tmp2[,c("True", "Prediction")] %>% table

# run Rosetta again to find further differences that might help to distinguish between those objects
ros_tmp <- rosetta(df_tmp1[,1:(dim(df_tmp1)[2] - 1)], discrete = T, underSample = T, reducer = "Johnson")
ros_tmp$quality

rec_tmp <- recalculateRules(df_tmp1[,1:(dim(df_tmp1)[2] - 1)], ros_tmp$main, discrete = T)
viewRules(head(rec_tmp, 20))
viewRules(rec_tmp[rec_tmp$decision == 1,])
viewRules(rec_tmp[rec_tmp$decision == 0,])


# Analysis of df_tmp2 -> check how many different values per feature
discern <- apply(df_tmp2[,1:(dim(df_tmp2)[2]-2)], 2, function(x){unique(x) %>% length}) %>% unlist
discern



