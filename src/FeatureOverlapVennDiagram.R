# This file contains steps to create a visualization to compare to rule tables
# created from R.Rosetta. The visualization should be a Venn Diagram that 
# shows the name of the common features and also their importance and the 
# associated decision

# Set working directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
set.seed(0)

# load necessary packages
library(dplyr)
library(VennDiagram)
library(ggplot2)
library(ggVennDiagram)

# given two Rosetta models (rule tables)
model1 <- rec[rec$pValue <= 0.05, ] # from RosettaOnNNLabels.R
model2 <- rec_rules[rec_rules$pValue <= 0.05, ] # from Rosetta.R

# extract a data frame for both rule table containing: feature, decision, p-Value
df1 <- model1[,c("features", "decision", "pValue")]
df2 <- model2[,c("features", "decision", "pValue")]

# split features in single rows (i.e., only one feature per row)
splitFeatures <- function(row){
  tmp <- strsplit(row[1], ',') %>% unlist %>% unname

  # no need for change if only one feature in row
  if (length(tmp) == 1){
    return(row)
  }
  # split rows if more than one feature per row
  df_tmp <- data.frame(matrix(NA, nrow = length(tmp), ncol = length(row)))
  df_tmp[,1] <- tmp
  df_tmp[,2:dim(df_tmp)[2]] <- matrix(row[2:length(row)], ncol = (dim(df_tmp)[2] -1), nrow = length(tmp), byrow = T)
  return(df_tmp %>% t)
}


df1 <- data.frame(matrix(apply(df1, 1, splitFeatures) %>% unlist %>% unname,
                         ncol = 3,
                         byrow = T))


df2 <- data.frame(matrix(apply(df2, 1, splitFeatures) %>% unlist %>% unname,
                         ncol = 3,
                         byrow = T))

# remove rows when feature-decision pair is duplicated
# i.e., keep most significant row if a duplication appear
df1 <- df1[!duplicated(df1[,1:2]),]
df2 <- df2[!duplicated(df2[,1:2]),]

# give meaningful column names
names(df1) <- c("Features", "Decision", "pValue")
names(df2) <- c("Features", "Decision", "pValue")

# create a Venn Diagram with the properties that are wished for
venn.diagram(x = list(df1$Features, df2$Features), 
             filename = "VennTest",
            col = c("darkgreen", "darkred"), 
            fill = c("lightgreen", "pink"), 
            alpha = c(0.5, 0.5),
            cat.col = c("darkgreen", "darkred"), 
            cat.cex = 1.5,
            #cex = c(df1$pValue, df2$pValue), 
            #area1 = nrow(df1),
            #area2 = nrow(df2),
            category.names = c("Model 1", "Model 2"))


# different approach for Venn Diagramm
y = Venn(list(df1$Features, df2$Features)) %>% process_data
ggplot() +
  geom_sf(aes(fill = count), data = venn_region(y)) +
  geom_sf(aes(color = id), data = venn_setedge(y), show.legend = FALSE) +
  geom_sf_text(aes(label = name), data = venn_setlabel(y)) +
  geom_sf_label(aes(label = count), data = venn_region(y)) +
  theme_void()
