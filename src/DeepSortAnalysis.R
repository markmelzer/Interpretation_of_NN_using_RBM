# Flush r-studio and set working directory 
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
set.seed(0)

# load necessary packages
library(dplyr)
library(caTools)
library(party)


# load results
res <- read.csv("../data/scDeepSort/human_Blood_human_Blood3223_data.csv")

# load truth
cell_types <- read.csv("../data/scDeepSort/human_Blood3223_celltype.csv")

# check that cells are ordered in same way
if(sum(cell_types$Cell == res$index) != dim(res)[1]){
  errorCondition("Cells do not match by index")
}

# assess quality of prediction
sum(res$cell_type == cell_types$Cell_type)/dim(cell_types)[1]
table(cell_types$Cell_type, res$cell_type)

# load gene expression data
dat <- read.csv("../data/scDeepSort/human_Blood3223_data.csv")
row.names(dat) <- dat[,1]
dat <- dat[,-1]

# check that cell index matches the ones from the predictions
if(sum(names(dat) == res$index) != dim(res)[1]){
  errorCondition("Cells do not match by index")
}

# transpose and add label at the end
df <- data.frame(t(dat))
df$Prediction <- factor(res$cell_type)

# split data into train and test set
sample_data <- sample.split(1:dim(df)[1], SplitRatio = 0.7)
train_data <- subset(df, sample_data == T)
test_data <- subset(df, sample_data == F)

# create control object for decision tree
ctl <- ctree_control(testtype = "Bonferroni", minsplit = 25, minbucket = 15, stump = F)

# set weights (increase for APC's to allow leave containing only APCs)
wgt <- rep.int(1, dim(train_data)[1])
wgt[train_data$Prediction == "Antigen presenting cell"] = 4

# create decision tree
model <- ctree(Prediction~., train_data, control = ctl, weights = wgt)
plot(model)

# predict test set
predict_model <- predict(model, test_data)
tab <- table(test_data$Prediction, predict_model)
sum(diag(tab))/sum(tab)
# 0.898


# find misclassified objects
misC <- (res$cell_type != cell_types$Cell_type) %>% which

table(res$cell_type[misC], cell_types$Cell_type[misC])

# cell identifiers of misclassified cells
cells <- cell_types$Cell[misC]

# their position in model
misC.in.tree <- model@responses@variables %>% row.names %in% cells %>% which

# check in which leave are unexpectedly many misclassifications

tab.leaves <- model@where %>% table # 2256 objects in leaves in total
tab.misc.leaves <- model@where[misC.in.tree] %>% table # 158 misclassifications

y = lapply(tab.leaves, function(x){158 * x / 2256})

# compare expected vs actual number
y %>% unlist
tab.misc.leaves

# find misclassified objects by misclassification

# confusion matrix
confusion_leaf <- function(x){
  return(table(train_data$Prediction[model@where == x], 
               cell_types$Cell_type[cell_types$Cell %in% row.names(train_data)[model@where == x]]))
}

confusion_leaf(7)
confusion_leaf(10)
confusion_leaf(14)
confusion_leaf(25)
confusion_leaf(27)
confusion_leaf(32)


# create data frame of objects in leafs
df_leaf <- function(x){
  tmp <- df[row.names(train_data)[model@where == x],  (1:dim(train_data)[2]-1)]
  tmp$Outcome <- factor(cell_types$Cell_type[cell_types$Cell %in% (row.names(train_data)[model@where == x])])
  return(tmp)
}


df_7  <- df_leaf(7)
df_10 <- df_leaf(10)
df_14 <- df_leaf(14)
df_25 <- df_leaf(25)
df_27 <- df_leaf(27)
df_32 <- df_leaf(32)

# create decision tree for these data frames
ctl2 <- ctree_control(testtype = "Bonferroni", minsplit = 3, minbucket = 2)

model.7 <- ctree(Outcome~., df_7, controls = ctl2)
model.10 <- ctree(Outcome~., df_10, controls = ctl2)
model.14 <- ctree(Outcome~., df_14, controls = ctl2)
model.25 <- ctree(Outcome~., df_25, controls = ctl2)
model.27 <- ctree(Outcome~., df_27, controls = ctl2)
model.32 <- ctree(Outcome~., df_32, controls = ctl2)


# plot decision trees
plot(model.7)
plot(model.10)
plot(model.14)
plot(model.25)
plot(model.27)
plot(model.32)

# find in which leafs the misclassified objects are now 
misC.in.tree7 <- (model.7@responses@variables %>% row.names) %in% cells %>% which
tab.misc.leaves.7 <- model.7@where[misC.in.tree7] %>% table
scLabel.misC.leaf7 <- model@responses@variables[(model.7@responses@variables %>% row.names)[(model.7@responses@variables %>% row.names %in% cells)],] %>% as.vector()
cbind(model.7@responses@variables[misC.in.tree7,] %>% as.vector, scLabel.misC.leaf7, model.7@where[misC.in.tree7])

misC.in.tree10 <- model.10@responses@variables %>% row.names %in% cells %>% which
tab.misc.leaves.10 <- model.10@where[misC.in.tree10] %>% table
scLabel.misC.leaf10 <- model@responses@variables[(model.10@responses@variables %>% row.names)[(model.10@responses@variables %>% row.names %in% cells)],] %>% as.vector()
cbind(model.10@responses@variables[misC.in.tree10,] %>% as.vector, scLabel.misC.leaf10, model.10@where[misC.in.tree10])

# trivial for the other trees

# ============= Check distribution of gene expression values ==================

# create density plots for the first three nodes and add line for split value

library(ggplot2)

df.FTL <- as.data.frame(cbind(cell_types$Cell_type, df[,"FTL"]))
names(df.FTL) <- c("CellType", "FTL")
df.FTL$FTL <- as.numeric(df.FTL$FTL)
ggplot(data = df.FTL, aes(x=CellType, y=FTL, fill = CellType), legend = F) +
  geom_violin() +
  geom_hline(yintercept = 3.756) +
  theme(legend.position="none") +
  ylab("FTL Expression")

df.CD20 <- as.data.frame(cbind(cell_types$Cell_type[df$FTL <= 3.756], df[df$FTL <= 3.756,"MS4A1"]))
names(df.CD20) <- c("CellType", "CD20")
df.CD20$CD20 <- as.numeric(df.CD20$CD20)
ggplot(data = df.CD20, aes(x=CellType, y=CD20, fill = CellType)) +
  geom_violin() +
  geom_hline(yintercept = 1.703) +
  theme(legend.position="none") +
  ylab("MS4A1 Expression")


df.TYROBP <- as.data.frame(cbind(cell_types$Cell_type[df$FTL > 3.756], df[df$FTL > 3.756,"TYROBP"]))
names(df.TYROBP) <- c("CellType", "TYROBP")
df.TYROBP$TYROBP <- as.numeric(df.TYROBP$TYROBP)
ggplot(data = df.TYROBP, aes(x=CellType, y=TYROBP, fill = CellType)) +
  geom_violin() +
  geom_hline(yintercept = 2.08) +
  theme(legend.position="none") +
  ylab("TYROBP Expression")





