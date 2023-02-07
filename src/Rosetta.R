# Flush r-studio and set working directory 
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
set.seed(0)

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
load_data <- function(file_name, header=T, sep=',', colClass="character", decisionNumeric=T) {
  dat = read.csv(file_name, header = header, sep = sep, colClasses = colClass)
  if(decisionNumeric){
    dat[,dim(dat)[2]] = as.numeric(dat[,dim(dat)[2]])
  }
  return(dat)
}


dat <- load_data("NS1.csv")

# MCFS
MCFS_features <- function(dat, decision="Pathogenicity~."){
  mcfs <- mcfs(Pathogenicity~., dat)
  attr <- append(mcfs$RI[1:mcfs$cutoff_value,]$attribute, names(dat[length(dat)]))
  return(attr)
}

"Read MCFS results from Zeeshans file: rmcfs gets an heap memory error as more
  memory is necessary than allowed by Java"
#mcfs <- read.csv("mcfs_result_Zeeshan.txt")

attr <- MCFS_features(dat)

# get Rosetta input data table
df = dat[attr]

# output from Uppmax
umax = readRDS("df_for_rosetta")
umax = umax[,names(umax)!="id"]

# undersampling


# run Rosetta

run_Rosetta <- function(df, roc = T, clroc = 1, discrete = T, undersample = T, underSampleNum = 0, method = "Johnson"){
  ros <- rosetta(df, roc = roc, clroc = clroc, discrete = discrete, underSample = undersample, underSampleNum = underSampleNum, reducer = method)
  return(ros)
}
ros = run_Rosetta(umax, underSampleNum = 10)
x = data.frame(cbind(umax[,1:100], umax[,101]))
names(x)[length(x)] = "Pathogenicity"
rosx = run_Rosetta(x, underSampleNum = 10)

# recalculateRules?

get_significant_rules <- function(rules, df, discrete = T, cutoff = 0.05){
  rules <- recalculateRules(df, rules, discrete = discrete)
  return(rules[rules$pValue <= 0.05,])
}

sig_rules <- get_significant_rules(ros$main, umax)


viewRules(head(sig_rules[sig_rules$decision=="1",]))

viewRules(head(sig_rules[sig_rules$decision=="0",]))


sig_rules = sig_rules[viewRules(sig_rules)$length<=1,]


# function to write features as list
manipulate_output <- function(sig_rules) {
  copy = sig_rules[,1:10]
  for (i in 1:length(copy$features)){
    new = c()
    old = unlist(strsplit(unlist(copy$features[i]), ','))
    lev = unlist(strsplit(unlist(copy$levels[i]), ','))
    for (j in old) {
      num = substr(j, 2, nchar(j))
      new = append(new, as.integer(num))
    }
    copy$features[i] = list(new)
    copy$levels[i] = list(lev)
  }
  return(copy)
}

copy = manipulate_output(sig_rules)


# function to write rules in hash map with position as key
get_hash_map <- function(sig_rules) {
  copy = sig_rules[,1:10]
  h <- hash()
  for (i in 1:length(copy$features)){
    old = unlist(strsplit(unlist(copy$features[i]), ','))
    for (j in old) {
      num = substr(j, 2, nchar(j))
      h[num] = append(h[[num]], i)
    }
  }
  return(h)
}

h = get_hash_map(sig_rules)



# function to get levels of position based on decision in rules
get_levels <- function(rules, hashmap, position, decision) {
  if (typeof(position) != "character"){
    position = as.character(position)
  }
  participating_rules <- rules[hashmap[[position]],]
  dec_rules <- participating_rules[participating_rules$decision == decision,]
  levels = c()
  position = as.integer(position)
  for (i in 1:length(dec_rules)){
    pos = which(unlist(dec_rules$features) == position)
    levels = append(levels, unlist(dec_rules$levels)[pos])
  }
  return(unique(levels))
}






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
"""
s = '['
for (i in quad) {
  s = paste(s, i)
  s = paste(s, ', ')
}
s = paste(s, ']')
"""
all = unique(append(single, append(double, triple)))

s = '['
for (i in sort(all)) {
  s = paste(s, i)
  s = paste(s, ', ')
}
s = paste(s, ']')
