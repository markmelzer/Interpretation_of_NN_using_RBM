# Flush r-studio and set working directory 
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
set.seed(0)

# load necessary packages
library(dplyr)
library(rmcfs)

# get data
dat <- read.csv("../data/Tubercolosis/Annotated_combined.csv", header = F, colClass = "character")

# remove all columns with only one amino acid -> cannot help classification
removed <- which((apply(dat, 2, function(x){length(table(x))}) == 1)) %>% unname
df <- dat[(apply(dat, 2, function(x){length(table(x))}) != 1) %>% unname]

# remove column containing mutation information
df <- df[,1:(dim(df)[2]-1)]

# change name of decision column
names(df)[dim(df)[2]] = "Resistance"

# run mcfs
res <- mcfs(Resistance~., df, cutoffPermutations = 50, projections = 50000)
res <- readRDS("../data/Tubercolosis/MCFS_output.RDS")
plot(res, type = "distances")

# extract attributes with RI greater than 0
attr <- res$RI$attribute[res$RI$RI > 0]

relev_df <- df[,c(attr, "Resistance")]


# save data
write.table(relev_df, "../data/Tubercolosis/Relevant_combi.csv", sep = ',', row.names = T, col.names = FALSE, quote = FALSE)
saveRDS(attr, "../data/Tubercolosis/Relevant_attributes.RDS")
saveRDS(res, "../data/Tubercolosis/MCFS_output.RDS")

