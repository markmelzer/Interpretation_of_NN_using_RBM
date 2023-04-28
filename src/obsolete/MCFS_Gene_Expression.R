# Flush r-studio and set working directory 
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
set.seed(0)


# load MCFS results from Uppmax

mcfs <- readRDS("../results/Exp3_GE/MCFS_output_GE.rds")
mcfs$cutoff_value
dim(mcfs$RI)


library(rmcfs)
plot(mcfs, type = "distances")
plot(mcfs, type = "ri")
plot(mcfs, type = 'features', size = 10)

# get data
dat <- as.data.frame(t(readRDS("../data/normalized_GE_data.RDS")))

# label data
load("../data/E-GEOD-68086-atlasExperimentSummary.Rdata")
dis <- experiment_summary@listData$rnaseq@colData$disease

subjects <- which(dis %in% c("normal", "breast carcinoma"))
#subjects <- which(dis %in% c("normal", "non-small cell lung carcinoma"))

label <- dis[subjects]
label[which(label != "normal")] = 1
label[which(label == "normal")] = 0
label = as.numeric(label)
dat$Outcome <- label

input <- dat[,c(mcfs$RI$attribute[1:100], "Outcome")]


# Rosetta

library(R.ROSETTA)

out <- rosetta(dt = input, discrete = F, discreteMethod = "EqualFrequency",
               discreteParam = 3, reducer = "Johnson", roc = T, clroc = 1,
               underSample = T)









