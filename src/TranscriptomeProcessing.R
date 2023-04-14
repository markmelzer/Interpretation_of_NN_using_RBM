# Flush r-studio and set working directory 
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
set.seed(0)

# get data
load("../data/E-GEOD-68086-atlasExperimentSummary.Rdata")

raw <- experiment_summary@listData$rnaseq@assays$data@listData$counts
batches <- experiment_summary@listData$rnaseq@colData$block

dim(raw)
# 65217 genes

# as in paper: exclude genes with less than five (non-normalized) read counts in all samples
raw <- raw[-which(apply(raw, 1, sum) < 5),]
# - 8419 genes

dim(raw)
# 56798 genes


# only use patients of breast cancer study
dis <- experiment_summary@listData$rnaseq@colData$disease

subjects <- which(dis %in% c("normal", "breast carcinoma"))
#subjects <- which(dis %in% c("normal", "non-small cell lung carcinoma"))

label <- dis[subjects]
label[which(label != "normal")] = 1
label[which(label == "normal")] = 0
label = as.numeric(label)
data <- raw[,subjects]

# get batches 
batch <- batches[subjects]



library(edgeR)
# weighted trimmed mean of M-values (TMM) --> library norm factor
y <- calcNormFactors(data)


# create design matrix
design <- model.matrix(~ batch + label)

# create DGE list from data
dge <- DGEList(counts = data, group = label, norm.factors = y)

# compute dispersion
disp <- estimateDisp(dge, design)
plotBCV(disp)

# =============================================================================
# fit generalized linear model
g <- glmFit(disp, design)
gt <- glmLRT(g, coef = length(colnames(design)))
t <- gt$table

# make color for plot
col = unname(abs(decideTests(gt)))
col[col == 0,] = "black"
col[col == 1,] = "red"

plot(2^t$logCPM, t$logFC, log = "x", col = col)


# =============================================================================
# fit (negative binomial) generalized linear model

glm <- glmQLFit(disp, design)
plotQLDisp(glm)

test <- glmQLFTest(glm, coef = 6)

tt <- test$table
#tt$decision <- unname(abs(decideTests(test)))
col = unname(abs(decideTests(test)))
col[col == 0,] = "black"
col[col == 1,] = "red"

plot(2^(tt$logCPM), tt$logFC, log = "x", col = col)


is.de <- decideTests(test, p.value=0.05)
summary(is.de)



# compare with EBI results
truth = read.csv("C:/Users/Mark/Downloads/E-GEOD-68086-analytics.tsv", sep = '\t', header = T)
FC = truth$X.breast.carcinoma..vs..normal..log2foldchange
res <- tt$logFC

qqplot(FC, res)
plot(2^(tt$logCPM), FC, log = "x", col = col)


# save fitted values to csv and RDS
write.table(test$fitted.values, file = "../data/normalized_GE_data.csv",
            quote = F, sep = ',', row.names = T, col.names = T)
saveRDS(test$fitted.values, "../data/normalized_GE_data.RDS")
