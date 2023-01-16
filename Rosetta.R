# Flush r-studio and set working directory 
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
set.seed(42)

# install and load necessary packages
#install.packages("digest")
library(digest)
#install.packages("rJava")
#install.packages("rmcfs")
#install.packages("devtools")
#library(devtools)
#install_github("komorowskilab/R.ROSETTA", force = T)

library(rJava)
library(rmcfs)
library(R.ROSETTA)

# load data
dat <- read.csv("NS1.csv", header = TRUE, sep= ',',
                colClasses = "character")
dat$Pathogenicity <- as.numeric(dat$Pathogenicity)

# MCFS
"Read MCFS results from Zeeshans file: rmcfs gets an heap memory error as more
  memory is necessary than allowed by Java"
mcfs <- read.csv("mcfs_result_Zeeshan.txt")
attr <- append(mcfs$attribute, names(dat[length(dat)]))
df = dat[attr]

# undersampling


# run Rosetta
ros_john <- rosetta(df, roc = TRUE, clroc = 1, discrete = TRUE,
                    underSample = TRUE) 

# recalculateRules?

rules_john <- ros_john$main

viewRules(head(rules_john[rules_john$decision=="1",]))

viewRules(head(rules_john[rules_john$decision=="0",]))


rules <- recalculateRules(df, rules_john, discrete = TRUE)
sig_rules <- rules[rules$pValue <= 0.05,]

head(viewRules(sig_rules))

# get array for Julia
feat = c()
lengths = c()
for (i in 1:length(sig_rules$features)) {
  feat[i] = (strsplit(sig_rules$features[i], ','))
  lengths[i] = length(unlist(feat[i]))
}

single = c()
double = c()
triple = c()
quad = c()
more = c()
for (i in feat){
  tmp = unlist(i)
  if (length(i) == 1) {
    single[length(single) + 1] = as.integer(substr(tmp, 2, nchar(tmp)))
  }
  else if (length(i) == 2){
    for (j in tmp){
      double[length(double) + 1] = as.integer(substr(j, 2, nchar(tmp)))
    }
  }
  else if (length(i) == 3){
    for (j in tmp){
      triple[length(triple) + 1] = as.integer(substr(j, 2, nchar(tmp)))
    }
  }
  else if (length(i) == 4){
    for (j in tmp){
      quad[length(quad) + 1] = as.integer(substr(j, 2, nchar(tmp)))
      }
    } else {
        for (j in tmp){
          more[length(more) + 1] = as.integer(substr(j, 2, nchar(tmp)))
        }
      }
}




