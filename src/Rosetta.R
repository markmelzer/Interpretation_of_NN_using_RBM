# Flush r-studio and set working directory 
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
set.seed(0)

# load necessary packages

library(rmcfs)
library(R.ROSETTA)

# load data
dat = read.csv("../data/HPAIV_train_set.csv", header=F, sep=',', colClass="character")
dat[,dim(dat)[2]] = as.numeric(dat[,dim(dat)[2]])

table(dat$V250)
# 1686 low, 512 high pathogenic

# MCFS

mcfs <- mcfs(V250~., dat, projections = 5000)

plot(mcfs, type = "distances")

attr <- append(mcfs$RI[1:mcfs$cutoff_value,]$attribute, names(dat[length(dat)]))

mcfs <- readRDS("../results/Exp1_HPAIV/MCFS_result.RDS")


"Read MCFS results from Zeeshans file: rmcfs gets an heap memory error as more
  memory is necessary than allowed by Java"
#mcfs <- read.csv("mcfs_result_Zeeshan.txt")


# get Rosetta input data table
df = dat[attr]

# output from Uppmax
umax = readRDS("df_for_rosetta")
umax = umax[,names(umax)!="id"]

# undersampling


# run Rosetta


ros = rosetta(df, roc = T, clroc = 1, discrete = T, underSample = T, underSampleNum = 5, reducer = "Johnson")
ros <- readRDS("../results/Exp1_HPAIV/Ros_result_Johnson.RDS")

x = data.frame(cbind(umax[,1:100], umax[,101]))
names(x)[length(x)] = "Pathogenicity"
rosx = run_Rosetta(x, underSampleNum = 10)

# get significant rules

rec_rules <- recalculateRules(dat, ros$main, discrete = T)
sig_rules <- rec_rules[rec_rules$pValue <= 0.05,]



viewRules(head(sig_rules[sig_rules$decision=="1",]))

viewRules(head(sig_rules[sig_rules$decision=="0",]))

filtered_rules <- sig_rules[sig_rules$coverageLHS >= 0.1,]
filtered_rules <- filtered_rules[filtered_rules$supportRHS >= 200,]


# get features
feat <- getFeatures(filtered_rules)

high_feat_1 <- feat$features$'1'[feat$frequencies$'1' >= 3]
high_feat_0 <- feat$features$'0'[feat$frequencies$'0' >= 3]

feat_0 <- feat$features$'0'
feat_1 <- feat$features$'1'

# after getting NN_attr below:
which(NN_attr %in% feat_1)
which(NN_attr %in% feat_0)

disc_features <- union(feat$features$'1', feat$features$'0')

get_num <- function(string){
  return(as.numeric(substr(string, 2, nchar(string))))
}

pos <- unlist(lapply(disc_features, get_num))
sort(pos)


# assess quality
test <- read.csv("../data/HPAIV_test_set.csv", header = F, colClasses = "character")
test[,dim(test)[2]] = as.numeric(test[,dim(test)[2]])
table(test[,250])

test_df <- test[attr]
write.table(test_df, "../data/HPAIV_test_set_imp_feat.csv", quote = F, col.names = F, row.names = F, sep = '\t')

pred <- predictClass(dt = test[,1:249], rules = sig_rules, discrete = T, validate = T, defClass = test[,250], 
                     normalize = T)

table(pred$out$currentClass, pred$out$predictedClass)
table(pred$out$currentClass)



# get NN features, sorted by weight small to high
w = c(74, 64, 53, 106, 50, 93, 108, 52, 16, 65, 32, 33, 61, 95, 99, 90, 73, 18, 66, 49, 78, 84, 48, 101,
      63, 39, 43, 40, 35, 45, 17, 55, 77, 70, 29, 103, 67, 23, 110, 51, 72, 15, 24, 26, 94, 27, 62, 68,
      22, 56, 42, 107, 6, 86, 36, 44, 104, 13, 102, 20, 30, 5, 112, 92, 79, 7, 28, 85, 105, 41, 25, 57,
      31, 38, 111, 37, 96, 75, 98, 10, 81, 46, 91, 21, 54, 87, 14, 59, 1, 89, 97, 19, 109, 71, 80, 100,
      60, 58, 8, 69, 83, 82, 47, 34, 88, 76, 9, 11, 12, 2, 3, 4)
w = w[length(w):1L]
NN_attr = attr[w]


