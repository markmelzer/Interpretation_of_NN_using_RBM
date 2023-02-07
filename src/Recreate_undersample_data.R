# Flush r-studio and set working directory 
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
set.seed(0)

feat = "P6,P7,P14,P18,P21,P22,P23,P24,P25,P26,P27,P28,P33,P42,P44,P48,P54,P55,P56,P59,P60,P63,P67,P70,P71,P73,P74,P81,P84,P96,P97,P100,P104,P105,P108,P111,P113,P118,P121,P122,P124,P126,P127,P128,P137,P147,P150,P155,P156,P163,P168,P173,P176,P180,P181,P202,P203,P205,P208,P209,P215,P216,P217,P218,P222,P223,P224,P228,P232,P236,P237,Pathogenicity"
feat = unlist(strsplit(feat, ','))


# load data
dat <- read.csv("NS1.csv", header = TRUE, sep= ',',
                colClasses = "character")
dat$Pathogenicity <- as.numeric(dat$Pathogenicity)
df = dat[,feat]

ros <- rosetta(df, roc = TRUE, clroc = 1, discrete = TRUE, underSample = TRUE, underSampleNum = 100)

rules <- recalculateRules(df, ros$main, discrete = TRUE)
sig_rules <- rules[rules$pValue <= 0.05,]

head(viewRules(sig_rules))


rules_80acc = sig_rules[sig_rules$accuracyRHS >= 0.8,]
total_LP = table(df$Pathogenicity)[1]
total_HP = table(df$Pathogenicity)[2]
csc_50 = c()
for (i in 1:dim(rules_80acc)[1]){
  if(rules_80acc[i,]$decision == "1"){
    csc = (rules_80acc[i,]$accuracyRHS * rules_80acc[i,]$supportLHS)/total_HP
  } else{
    csc = (rules_80acc[i,]$accuracyRHS * rules_80acc[i,]$supportLHS)/total_LP
  }
  if(csc >= 0.5){
    csc_50 = append(csc_50, TRUE)
  }else{
    csc_50 = append(csc_50, FALSE)
  }
}
strong_rules = rules_80acc[csc_50,]
