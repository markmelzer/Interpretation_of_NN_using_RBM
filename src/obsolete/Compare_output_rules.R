zeeshan = read.csv("rosetta_result_Zeeshan.txt", sep = '\t', header = FALSE)


rules_80acc = sig_rules_u[sig_rules_u$accuracyRHS >= 0.8,]
total_LP = table(umax$Pathogenicity)[1]
total_HP = table(umax$Pathogenicity)[2]
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

r = strong_rules

rule = c()
for (i in 1:dim(r)[1]) {
  feat = unlist(strsplit(r[i,]$features, ','))
  lev = unlist(strsplit(r[i,]$levels, ','))
  tmp = ""
  for (j in 1:length(feat)){
    tmp = paste0(tmp, feat[j])
    tmp = paste0(tmp, '=')
    tmp = paste0(tmp, lev[j])
    tmp = paste0(tmp, ',')
  }
  rule[length(rule)+1] = tmp
}
mark = cbind(rule, r$decision, r$accuracyRHS, r$supportRHS)
mark = as.data.frame(mark)

?intersect
intersect(mark$rule, zeeshan$V1)
length(union(mark$rule, zeeshan$V1)) - length(intersect(mark$rule, zeeshan$V1))

# add quality processing of Zeeshan in paper