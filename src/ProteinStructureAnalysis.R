# Script to create a heatmap that compares structural elements of a protein with
# gaps in the consensus sequence of an MSA, and the rules of a Rosetta model

# Flush r-studio and set working directory 
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
set.seed(0)

# load necessary library
library(bio3d)
library(dplyr)
library(tidyverse)
library(Biostrings)

# load structure information for NS1
struct <- read.pdb("../data/NS1/3f5t.pdb")

# load aligned sequence (from AlignAASeq2MSA.R)
seq <- read.csv("../data/NS1/alignedSeq.fa")$x

# get gap positions
gap <- lapply(seq, function(x){x == "-"}) %>% unlist

# compute position of aligned sequence residues in original sequence
gap_count <- array(0, 249)
for(i in 1:249){
  gap_count[i] <- sum(gap[1:i])
}
true_pos <- c(1:249) - gap_count

# mark gaps by using 0
true_pos[gap]<- 0

# make array to save helix and sheet position
helix <- array(0, 249)

for (i in 1:length(struct$helix$start)) {
  start <- unname(struct$helix$start)[i]
  start <- start + gap_count[start]
  end <- unname(struct$helix$end)[i]
  end <- end + gap_count[end]
  helix[start:end] = array(1, (end-start+1))
}

sheet <- array(0, 249)

for (i in 1:length(struct$sheet$start)) {
  start <- unname(struct$sheet$start)[i]
  start <- start + gap_count[start]
  end <- unname(struct$sheet$end)[i]
  end <- end + gap_count[end]
  sheet[start:end] = array(1, (end-start+1))
}

# read the MSA
df <- read.csv("../data/NS1/NS1_H5_H7.csv", header = T, colClasses = "character")

# extract the columns of serotype H5N1
serotype <- unlist(lapply(df$id, function(x){unlist(strsplit(x, '[*]'))[5]}))
h5n1 <- df[which(serotype=="H5N1"),]
msa <- df[,2:(dim(df)[2]-1)]

# compute consensus sequence for positive and negative sequence
pos <- msa[df$Pathogenicity == "1",]
neg <- msa[df$Pathogenicity == "0",]

# change gap symbol to -
pos[pos == "?"] = "-"
neg[neg == "?"] = "-"

# load alignments
alignment_pos <- AAStringSet(apply(pos, 1, function(x){
  gsub('[, ]', '', toString(x))}) %>% unlist)

alignment_neg <- AAStringSet(apply(neg, 1, function(x){
  gsub('[, ]', '', toString(x))}) %>% unlist)

# compute consensus sequence
cons_pos <- consensusString(alignment_pos) %>% strsplit('') %>% unlist
cons_neg <- consensusString(alignment_neg) %>% strsplit('') %>% unlist

# get positions with consensus gaps
gaps_pos <- (cons_pos == "-") %>% as.numeric
gaps_neg <- (cons_neg == "-") %>% as.numeric

# get positions that are important according to R.ROSETTA
rec <- readRDS("../results/Exp1_HPAIV/FullRuleTable.RDS")

# extract features (i.e., positions) used in Rosetta Rules
pos_rules <- rec[rec$decision == 1,]
neg_rules <- rec[rec$decision == 0,]

pos_features <- lapply(pos_rules$features, function(x){strsplit(x, ",") %>% unlist}) %>% unlist %>% unique
neg_features <- lapply(neg_rules$features, function(x){strsplit(x, ",") %>% unlist}) %>% unlist %>% unique

pos_feat <- lapply(pos_features, function(x){(strsplit(x, "V") %>% unlist)[2]}) %>% unlist %>% as.numeric
neg_feat <- lapply(neg_features, function(x){(strsplit(x, "V") %>% unlist)[2]}) %>% unlist %>% as.numeric

both <- intersect(pos_feat, neg_feat)
only_pos <- setdiff(pos_feat, both)
only_neg <- setdiff(neg_feat, both)

# Heat map of consensus gaps and positions

dat2 <- data.frame(2*helix, 2*sheet, gaps_neg, gaps_pos)
names(dat2) <- c("Helices", "Sheets", "Consensus gaps LP", "Consensus gaps HP")

dt <- dat2 %>% 
  rownames_to_column() %>%
  gather(colname, value, -rowname)

dt$rowname <- as.numeric(dt$rowname)

ggplot(dt, aes(x = rowname, y = colname, fill = value)) +
  geom_tile(show.legend = F) +
  geom_vline(xintercept = only_neg, color = "red") +
  geom_vline(xintercept = only_pos, color = "green") +
  geom_vline(xintercept = both, color = "yellow") +
  scale_fill_gradientn(colors = c("gray", "black")) +
  ylab("") +
  xlab("MSA Position")

