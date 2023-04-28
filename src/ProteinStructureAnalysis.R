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
df <- read.csv("../data/NS1/NS1.csv", header = T, colClasses = "character")

# extract the columns of serotype H5N1
serotype <- unlist(lapply(df$id, function(x){unlist(strsplit(x, '[*]'))[5]}))
h5n1 <- df[which(serotype=="H5N1"),]
msa <- df[,2:(dim(df)[2]-1)]

# compute consensus sequence for positive and negative sequence
pos <- msa[df$Pathogenicity == "1",]
neg <- msa[df$Pathogenicity == "0",]

pos[pos == "?"] = "-"
neg[neg == "?"] = "-"

alignment_pos <- AAStringSet(apply(pos, 1, function(x){
  gsub('[, ]', '', toString(x))}) %>% unlist)

alignment_neg <- AAStringSet(apply(neg, 1, function(x){
  gsub('[, ]', '', toString(x))}) %>% unlist)

cons_pos <- consensusString(alignment_pos) %>% strsplit('') %>% unlist
cons_neg <- consensusString(alignment_neg) %>% strsplit('') %>% unlist

gaps_pos <- (cons_pos == "-") %>% as.numeric
gaps_neg <- (cons_neg == "-") %>% as.numeric

# Heat map of consensus gaps and positions

dat2 <- data.frame(2*helix, 2*sheet, gaps_neg, gaps_pos)
names(dat2) <- c("Helices", "Sheets", "Consensus gaps LP", "Consensus gaps HP")

dt <- dat2 %>% 
  rownames_to_column() %>%
  gather(colname, value, -rowname)

dt$rowname <- as.numeric(dt$rowname)

ggplot(dt, aes(x = rowname, y = colname, fill = value)) +
  geom_tile(show.legend = F) +
  geom_vline(xintercept = c(86), color = "red") +
  geom_vline(xintercept = c(55, 27, 70, 60, 18, 84, 14, 6, 7, 26), color = "green") +
  geom_vline(xintercept = c(85), color = "yellow") +
  scale_fill_gradientn(colors = c("gray", "black")) +
  ylab("") +
  xlab("MSA Position")



# vertical lines at positions used in important rules by Rosetta model