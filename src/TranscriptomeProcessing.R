# Flush r-studio and set working directory 
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
set.seed(0)

# get data
load("C:/Users/Mark/Downloads/E-GEOD-68086-atlasExperimentSummary.Rdata")

#raw <- read.csv("../data/GeneExpression/E-GEOD-68086-raw-counts.tsv", sep = '\t')
raw <- experiment_summary@listData$rnaseq@assays$data@listData$counts

dim(raw)
# 65217 genes

# as in paper: exclude gnes with less than five (non-normalized) read counts in all samples
raw <- raw[-which(apply(raw, 1, sum) <= 5),]
# - 8419 genes

dim(raw)
# 56798 genes


# only use patients of breast cancer study
breast = raw[,c(3:26,252:266)]
reference = raw[,c(115:170, 249:251)]
reference = reference[,-32]
dim(breast)
dim(reference)


data <- data.frame(cbind(breast, reference))
#data <- raw


# weighted trimmed mean of M-values (TMM) --> library norm factor
y <- calcNormFactors(data)

# labels: 1 for breast cancer, 0 for reference
label <- append(rep(1, dim(breast)[2]), rep(0, dim(reference)[2]))

# create DGE list from data
dge <- DGEList(counts = data, group = label, norm.factors = y)


# compute dispersion
disp <- estimateDisp(dge)



# fit (negative binomial) generalized linear model
glm <- glmQLFit(disp)

test <- glmQLFTest(glm)


tt <- topTags(test, n=Inf)$table
tt$decision <- unname(abs(decideTests(test)))

plot(2^tt$logCPM, tt$logFC, log = "x", col = unname(abs(decideTests(test))))

plot(2^tt$logCPM, tt$logFC, log = "x")


#================= test ===================================

high_exp_genes <- data.frame(test$fitted.values[rownames(tt[tt$logCPM >= 3,]),])


#tmp = t(rbind(high_exp_genes, label))
#write.csv(tmp, "high_exp_genes.csv")

# clustering: Pearson distances, Ward clustering

library(ClassDiscovery)
d <- distanceMatrix(high_exp_genes, "pearson")

cl <- hclust(d, method = "ward.D")


b = data.frame(test$fitted.values[,1:39])
c = data.frame(test$fitted.values[,40:97])

b_means = apply(b, 1, mean)
c_means = apply(c, 1, mean)

b_means[b_means == 0] = 10^-8
c_means[c_means == 0] = 10^-8


M = log2(b_means/c_means)
A = 1/2 * (b_means+c_means)


plot(A, M, log = "x", ylim = c(-10, 10), xlim = c(0.001, 10^7))

df = data.frame(M, 2^A)

ggplot(df, aes(x=A,y=M)) + 
  geom_point() +
  scale_x_log10() +
  ylim(-10, 10)











