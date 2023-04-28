# Flush r-studio and set working directory 
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
set.seed(0)

# load necessary packages
library(dplyr)
library(R.ROSETTA)

# load data
dat <- read.csv("../AADataNewLabel.csv", header = F, colClasses = "character")

# extract true labels and remove from data frame
truth <- dat[,(dim(dat)[2] - 1)]
df <- dat[,c(1:(dim(dat)[2] - 2), dim(dat)[2])]

# run Rosetta on NN output
ros <- rosetta(df, discrete = T, underSample = T)
ros$quality
rec <- recalculateRules(df, ros$main, discrete = T)


# extract wrongly classified objects
wrongObj <- which(df[,dim(df)[2]] != truth)

# get supportSet as list
supportSet <- lapply(rec$supportSetRHS, function(x){as.numeric(unlist(strsplit(x, ",")))})

rulesOfSig <- unlist(lapply(supportSet, function(x){any(wrongObj %in% x)}))

supportSet[rulesOfSig]


# get an overview over the wrongly classified objects

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

for(i in 1:length(wrongObj)){
  print(paste0("Object Nr ", wrongObj[i]))
  tmp = which(unlist(lapply(supportSet, function(x){(wrongObj[i] %in% x)})))
  print(paste0("Appears in ", length(tmp), " rules"))
  print(df[i, dim(df)[2]])
  print(viewRules(rec[tmp,]))
  f = getFeatures(rec[tmp,])
  #print("Containg those features (0/1):")
  #print(f$features$'0')
  #print(f$features$'1')
  print("++++++++++++++++++++++++++++", quote = F)
  
}

# investigate the 6 wrong objects appearing in many rules
obj <- c(68, 125, 179, 317, 423, 488)
uni <- c()
inter <- c()
for(i in obj){
  tmp = which(unlist(lapply(supportSet, function(x){(i %in% x)})))
  f = getFeatures(rec[tmp,])
  print(f)
  print(union(f$features$'0', f$features$'1'))
  if(length(inter) == 0){
    # rule contains either decision 0 or 1
    inter <- union(f$features$'0', f$features$'1')
  } else{
    inter <- intersect(inter, union(f$features$'0', f$features$'1'))
  }
  uni <- union(uni, union(f$features$'0', f$features$'1'))
}

setdiff(uni, inter)
View(dat[obj, c(inter, "V250", "V251")])
View(dat[obj, c(setdiff(uni, inter), "V250", "V251")])

apply(dat[,inter], 2, table)

sum(unlist(lapply(supportSet, function(x){any(obj %in% x)})))

# check if the objects are in the LHS support set of other rules

LHSSupportSet <- lapply(rec$supportSetLHS, function(x){as.numeric(unlist(strsplit(x, ",")))})

LHSrulesOfSig <- unlist(lapply(LHSSupportSet, function(x){any(obj %in% x)}))
RHSrulesOfSig <- unlist(lapply(supportSet, function(x){any(obj %in% x)}))

# check all objects that have the same values in the intersection features (inter)
df_tmp <- filter(dat, V100 == "I" & V104 == "S" & V111 == "E" & V121 == "I"
                 & V14 == "Y" & V168 == "A" & V18 == "I" & V22 == "L"
                 & V222 == "G" & V236 == "R" & V6 == "I" 
                 & (V67 == "N" | V67 == "D") & V7 == "T" & V70 == "K"
                 & V85 == "A")

# run Rosetta on this new df
tmp_ros <- rosetta(df_tmp[,1:250], underSample = T, discrete = T, reducer = "Genetic")
tmp_ros$quality
tmp_rec <- recalculateRules(df_tmp[,1:250], tmp_ros$main, discrete = T)

# investigate the other 10 wrongly classified objects
wrong <- setdiff(wrongObj, obj)
wrong_fp <- wrong[which(dat[wrong, 250] == 0)]
wrong_fn <- setdiff(wrong, wrong_fp)

tt = wrong_fp
for(i in 1:length(tt)){
  print(paste0("Object Nr ", tt[i]))
  tmp = which(unlist(lapply(supportSet, function(x){(tt[i] %in% x)})))
  print(paste0("Appears in ", length(tmp), " rules"))
  #print(df[i, dim(df)[2]])
  print(viewRules(rec[tmp,]))
  f = getFeatures(rec[tmp,])
  print("Containg those features (0/1):")
  print(f$features$'0')
  print(f$features$'1')
  print("++++++++++++++++++++++++++++", quote = F)
}
viewRules(rec[which(unlist(lapply(LHSSupportSet, function(x){any(321 %in% x)}))),]) # check for different objects in wrong
df_tmp2 <- filter(dat, V221 == "R" & V146 == "I" & V44 == "G")

#======================== Rosetta on Train1 Data ==============================
train = read.csv("../data/HPAIV_Train1.csv", header=F, sep=',', colClass="character")
train[,dim(train)[2]] = as.numeric(train[,dim(train)[2]])
ros_train = rosetta(train, roc = T, clroc = 1, discrete = T, underSample = T, reducer = "Johnson")
rec_train = recalculateRules(train, ros_train$main)


train_feat <- getFeatures(ros_train$main[ros_train$main$pValue <= 0.05, ])
nn_feat <- getFeatures(ros$main[ros$main$pValue <= 0.05, ])

intersect(train_feat$features$'0', nn_feat$features$'0') %>% length
intersect(train_feat$features$'1', nn_feat$features$'1') %>% length

setdiff(nn_feat$features$'0', train_feat$features$'0') %>% sort
setdiff(nn_feat$features$'1', train_feat$features$'1') %>% sort

nn_feat$features$'0' %>% length

intersect(intersect(train_feat$features$'0', nn_feat$features$'0'), intersect(train_feat$features$'1', nn_feat$features$'1'))


# ======================= use findings on Test Set ============================
test <- read.csv("../data/HPAIV_Test.csv", header=F, sep=',', colClass="character")

err1 <- filter(test, V89 == "V")
err2 <- filter(test, V100 == "I" & V104 == "S" & V111 == "E" & V121 == "I"
              & V14 == "Y" & V168 == "A" & V18 == "I" & V22 == "L"
              & V222 == "G" & V236 == "R" & V6 == "I" 
              & (V67 == "N" | V67 == "D") & V7 == "T" & V70 == "K"
              & V85 == "A")

erridx1 <- which(test$V89 == "V")
experr1 <- 6 * (which(test$V89 == "V" & test$V250 == 0) %>% length)/(which(dat$V89 == "V" & dat$V251 == 0) %>% length)

erridx2 <- which(test$V100 == "I" & test$V104 == "S" & test$V111 == "E" & test$V121 == "I" &
        test$V14 == "Y" & test$V168 == "A" & test$V18 == "I" & test$V22 == "L" &
        test$V222 == "G" & test$V236 == "R" & test$V6 == "I" &
        (test$V67 == "N" | test$V67 == "D") & test$V7 == "T" &
        test$V70 == "K" & test$V85 == "A") # expect 6/74 * 44 wrongs (3.57) ???
experr2 <- 6 * (erridx2 %>% length)/sum(df_tmp$V250 == 0)

# run classification in Julia to check

dat2 <- read.csv("../AATestDataNewLabel.csv", header = F, colClasses = "character")

# extract true labels and remove from data frame
truth2 <- dat2[,(dim(dat2)[2] - 1)]
df2 <- dat2[,c(1:(dim(dat2)[2] - 2), dim(dat2)[2])]

# get wrong objects

wrongObj2 <- which(truth2 != df2[,dim(df2)[2]])

wrongObj2 %in% erridx1 %>% sum # 11 (exp: 4.12)
wrongObj2 %in% erridx2 %>% sum #  2 (exp: 3.88)

intersect(wrongObj2[wrongObj2 %in% erridx1], wrongObj2[wrongObj2 %in% erridx2]) # empty intersection as in Train2 before
wrongRest <- setdiff(wrongObj2, union(wrongObj2[wrongObj2 %in% erridx1], wrongObj2[wrongObj2 %in% erridx2]))

# run Rosetta on NN test output
library(R.ROSETTA)

ros2 <- rosetta(df2, discrete = T, underSample = T)
ros2$quality

rec2 <- recalculateRules(df2, ros2$main, discrete = T)

# check out the remaining 4 wrong objects
supportSet2 <- lapply(rec2$supportSetRHS, function(x){as.numeric(unlist(strsplit(x, ",")))})
rulesOfSig2 <- unlist(lapply(supportSet2, function(x){any(wrongRest %in% x)}))
unlist(lapply(supportSet2, function(x){sum(wrongRest %in% x)}))

dat2[wrongRest, 250:251]


