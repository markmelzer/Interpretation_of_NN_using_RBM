# Flush r-studio and set working directory 
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
set.seed(0)

# load necessary packages
library(R.ROSETTA)

# load data
dat = read.csv("../data/HPAIV_Train1.csv", header=F, sep=',', colClass="character")
dat[,dim(dat)[2]] = as.numeric(dat[,dim(dat)[2]])

table(dat$V250)
# 1153 low, 347 high pathogenic

# run Rosetta
ros_train = rosetta(dat, roc = T, clroc = 1, discrete = T, underSample = T, reducer = "Johnson")
#ros <- readRDS("../results/Exp1_HPAIV/Ros_result_Johnson.RDS")

# get significant rules
rec_rules <- recalculateRules(df, ros_train$main, discrete = T)
sig_rules <- rec_rules[rec_rules$pValue <= 0.05,]


# show rules
viewRules(head(sig_rules[sig_rules$decision=="1",]))
viewRules(head(sig_rules[sig_rules$decision=="0",]))


# get features
feat <- getFeatures(filtered_rules)

high_feat_1 <- feat$features$'1'[feat$frequencies$'1' >= 3]
high_feat_0 <- feat$features$'0'[feat$frequencies$'0' >= 3]

feat_0 <- feat$features$'0'
feat_1 <- feat$features$'1'


# assess quality on test set
test <- read.csv("../data/HPAIV_test_set.csv", header = F, colClasses = "character")
test[,dim(test)[2]] = as.numeric(test[,dim(test)[2]])
table(test[,250])

test_df <- test[attr]
write.table(test_df, "../data/HPAIV_test_set_imp_feat.csv", quote = F, col.names = F, row.names = F, sep = '\t')

pred <- predictClass(dt = test[,1:249], rules = sig_rules, discrete = T, validate = T, defClass = test[,250], 
                     normalize = T)

table(pred$out$currentClass, pred$out$predictedClass)
table(pred$out$currentClass)

