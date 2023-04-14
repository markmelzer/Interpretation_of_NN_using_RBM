# Flush r-studio and set working directory 
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
set.seed(0)

dat <- read.csv("../AADataNewLabel.csv", header = F, colClasses = "character")

# extract true labels and remove from data frame
truth <- dat[,(dim(dat)[2] - 1)]
df <- dat[,c(1:(dim(dat)[2] - 2), dim(dat)[2])]

# run Rosetta on NN output
library(R.ROSETTA)

ros <- rosetta(df, discrete = T)
ros$quality
rec <- recalculateRules(df, ros$main, discrete = T)

rec$supportSetRHS

# extract wrongly classified objects
wrongObj <- which(df[,dim(df)[2]] != truth)

# get supportSet as list
supportSet <- lapply(rec$supportSetRHS, function(x){as.numeric(unlist(strsplit(x, ",")))})

rulesOfSig <- unlist(lapply(supportSet, function(x){any(wrongObj %in% x)}))

supportSet[rulesOfSig]


for(i in 1:length(supportSet)){
  if(rulesOfSig[i]){
    x = unlist(supportSet[i])
    print(viewRules(rec[i,]))
    print("Proportion of support objects in wrong objects")
    print(sum(x %in% wrongObj)/(length(wrongObj)))
    print("Number and proportion of wrong objects in support set")
    print(sum(wrongObj %in% x))
    print(sum(wrongObj %in% x)/(length(x)))
  }
}

viewRules(rec[rulesOfSig,])
rec[rulesOfSig,]
