# Flush r-studio and set working directory 
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
set.seed(0)

# allocate java memory
options(java.parameters = "-Xmx128g")

# get data
dat <- t(readRDS("../data/normalized_GE_data.RDS"))

# get labels for patients
load("../data/E-GEOD-68086-atlasExperimentSummary.Rdata")
dis <- experiment_summary@listData$rnaseq@colData$disease

subjects <- which(dis %in% c("normal", "breast carcinoma"))

label <- dis[subjects]
label[which(label != "normal")] = 1
label[which(label == "normal")] = 0
label = as.numeric(label)

# check that names of patients are the same
rownames(experiment_summary@listData$rnaseq@colData)[subjects] == row.names(dat)

# remove unneccesary variables
rm(subjects)
rm(dis)
rm(experiment_summary)

# add label to data matrix
dat <- as.data.frame(dat)
dat$Cancer <- label

# MCFS
library(rmcfs)
res <- mcfs(Cancer ~ ., dat)






