# Flush r-studio and set working directory 
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
set.seed(0)


# load necessary packages
library(R.ROSETTA)
library(arules)


# load necessary data
test_MSA <- read.csv("../data/HPAIV_test_set.csv", sep=',', colClass="character", header = F)
test_MSA[,dim(test_MSA)[2]] = as.numeric(test_MSA[,dim(test_MSA)[2]])


# load values of NN nodes
node_val <- read.csv("../results/output_node_values.csv", header = F)
node_val <- node_val[2:dim(node_val)[1],]
node_val$V1 <- as.numeric(node_val$V1)

# to test: use first layer (but not input)

# extract neuron values for first layer
positions1 = 2+(0:233)*6
positions2 = 3+(0:233)*6
positions3 = 4+(0:233)*6
positions4 = 5+(0:233)*6




layer1 = node_val[positions1,]
layer2 = node_val[positions2,]
layer3 = node_val[positions3,]
layer4 = node_val[positions4,]

# remove columns with less values than cuts
f <- function(x){
  return(length(unique(x))!=1)
}

layer3 <- layer3[apply(layer3, 2, f)]

# discretize
discretize_data <- function(data, params = list(method = "interval", breaks = 2, 
                                                labels = c("1", "2"))){
  return(discretizeDF(data, default = params))
}
df1 = discretize_data(layer1)
df2 = discretize_data(layer2)
df3 = discretize_data(layer3, params = list(method = "cluster", breaks = 2, 
                                            labels = c("1", "2")))
df4 = discretize_data(layer4)




# add decision at the end
df1 = data.frame(cbind(df1, test_MSA[,dim(test_MSA)[2]]))
nam <- c(names(df1[1:length(names(df1))-1]), "Pathogenicity")
names(df1) <- nam
rm(nam)

df2 = data.frame(cbind(df2, test_MSA[,dim(test_MSA)[2]]))
nam <- c(names(df2[1:length(names(df2))-1]), "Pathogenicity")
names(df2) <- nam
rm(nam)

df3 = data.frame(cbind(df3, test_MSA[,dim(test_MSA)[2]]))
nam <- c(names(df3[1:length(names(df3))-1]), "Pathogenicity")
names(df3) <- nam
rm(nam)

df4 = data.frame(cbind(df4, test_MSA[,dim(test_MSA)[2]]))
nam <- c(names(df4[1:length(names(df4))-1]), "Pathogenicity")
names(df4) <- nam
rm(nam)



# run Rosetta
ros0 <- rosetta(test_MSA, clroc = 1, discrete = T, reducer = "Johnson")
ros1 <- rosetta(df1, clroc = 1, discrete = T, reducer = "Johnson")
ros2 <- rosetta(df2, clroc = 1, discrete = T, reducer = "Johnson")
ros3 <- rosetta(df3, clroc = 1, discrete = T, reducer = "Johnson")
ros4 <- rosetta(df4, clroc = 1, discrete = T, reducer = "Johnson")


ros0$quality
ros1$quality
ros2$quality
ros3$quality
ros4$quality

