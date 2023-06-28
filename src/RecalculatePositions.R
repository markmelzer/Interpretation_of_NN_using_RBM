# Flush r-studio and set working directory 
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
set.seed(0)

# load necessary packages
library(dplyr)

# ====== Recalculate positions for MTB ========================================

# get data
dat <- read.csv("../data/Tubercolosis/Annotated_combined.csv", header = F, colClass = "character")

# recalculation of positions
katG <- dat[,271:dim(dat)[2]]

# for position 191
gap191 <- apply(katG[,1:190], 1, function(x){table(x %>% unlist %>% unname)["-"]})
# suggest shift by 3 positions, check for gap in these positions
gap194 <- apply(katG[,191:194], 1, function(x){table(x %>% unlist %>% unname)["-"]})
# thus, position 191 of katG is shifted to 194 in the MSA, V464 of the concatenated MSA

# for position 315
gap315 <- apply(katG[,1:314], 1, function(x){table(x %>% unlist %>% unname)["-"]})

# suggest for shift of 9 (in three cases 10 positions)
gap325 <- apply(katG[,315:325], 1, function(x){table(x %>% unlist %>% unname)["-"]})
# thus, P315 shifted by 9 in MSA to 324 (concatenated to 594)
# in 5 cases by 10 positions, check which objects
out1 <- which(gap325 == 1)
out2 <- which(gap315 != 9)

# 15, 23, 83, 351, 610, 951 (check if objects are in Train2 and potentially misclassified by NN)
obj_id <- read.csv("../data/Tubercolosis/Tb_Train2_objects.csv", header = F)
out1 %in% obj_id$V1
out2 %in% obj_id$V1

# ====== Recalculate positions for MTB ========================================

# gaps for NS1 AIV
dat <- read.csv("../data/NS1/NS1_H5_H7.csv", header = T, colClasses = "character")
df <- dat[,2:(dim(dat)[2]-1)]
hp <- df[dat$Pathogenicity == 1,]
lp <- df[dat$Pathogenicity == 0,]

gap <- apply(df, 1, function(x){table(x %>% unlist %>% unname)["?"]})

# check how many gaps in front of position 84

hp.gap.84 <-apply(hp[,1:84], 1, function(x){table(x %>% unlist %>% unname)["?"]})
sum(hp.gap.84 != 5)
hp.gap.84 %>% table

lp.gap.84 <-apply(lp[,1:84], 1, function(x){table(x %>% unlist %>% unname)["?"]})
sum(lp.gap.84 != 5)
lp.gap.84 %>% table

# check how many gaps in front of position 228

hp.gap.228 <-apply(hp[,1:228], 1, function(x){table(x %>% unlist %>% unname)["?"]})
hp.gap.228 %>% table

lp.gap.228 <-apply(lp[,1:228], 1, function(x){table(x %>% unlist %>% unname)["?"]})
lp.gap.228 %>% table

# check differences between HP and LP sequences

hp.gap <-apply(hp, 1, function(x){table(x %>% unlist %>% unname)["?"]})
hp.gap %>% table

lp.gap <-apply(lp, 1, function(x){table(x %>% unlist %>% unname)["?"]})
lp.gap %>% table


