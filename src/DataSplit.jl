# split data into two training and one test set
using DelimitedFiles
using Random

# NS1 (MTB below)

# set seed to reproduce split
Random.seed!(0)

# read file and check dimensions
file = readdlm("../data/NS1/NS1_H5_H7.csv", ',')
size(file)

# remove header and id --> unnecessary for training
file = file[2:end,2:end]

# randomly reorder objects 
file= file[shuffle(1:end), :]

# check class distribution for imbalance
sum(file[2:end, end])
size(file)[1] - sum(file[2:end, end]) - 1 # -1 for header 

# split: NN probably needs more training data than RBM
# NS1: 549-250-101

# get positive and negative objects separately
pos = file[findall(isone, file[2:end, end]).+1,:]
neg = file[findall(iszero, file[2:end, end]).+1,:]

# split positive and negative objects

pos_a = pos[1:331,:]
pos_b = pos[332:482, :]
pos_c = pos[483:end, :]

neg_a = neg[1:217,:]
neg_b = neg[218:316, :]
neg_c = neg[317:end, :]

# combine positive and negative objects for the three data sets
train1 = vcat(pos_a, neg_a)[shuffle(1:end), :]
train2 = vcat(pos_b, neg_b)[shuffle(1:end), :]
test = vcat(pos_c, neg_c)[shuffle(1:end), :]

# write to files
writedlm("./data/NS1/NS1_H5_H7_Train1.csv", train1, ',')
writedlm("./data/NS1/NS1_H5_H7_Train2.csv", train2, ',')
writedlm("./data/NS1/NS1_H5_H7_Test.csv", test, ',')



# MTB

# set seed to reproduce split
Random.seed!(0)

# read file and check dimensions
file = readdlm("../data/Tubercolosis/Annotated_combined.csv", ',')
size(file)


# randomly reorder objects 
file= file[shuffle(1:end), :]

# check class distribution for imbalance
sum(file[:, end-1])
size(file)[1] - sum(file[:, end-1]) # -1 for header 

# split: NN probably needs more training data than RBM
# NS1: 632-209-148

# get positive and negative objects separately
pos = file[findall(isone, file[:, end-1]),:]
neg = file[findall(iszero, file[:, end-1]),:]

# split positive and negative objects

pos_a = pos[1:85,:]
pos_b = pos[86:125, :]
pos_c = pos[126:end, :]

neg_a = neg[1:547,:]
neg_b = neg[548:716, :]
neg_c = neg[617:end, :]

# combine positive and negative objects for the three data sets
train1 = vcat(pos_a, neg_a)[shuffle(1:end), :]
train2 = vcat(pos_b, neg_b)[shuffle(1:end), :]
test = vcat(pos_c, neg_c)[shuffle(1:end), :]

# write to files
writedlm("./data/Tubercolosis/Tb_Train1.csv", train1, ',')
writedlm("./data/Tubercolosis/Tb_Train2.csv", train2, ',')
writedlm("./data/Tubercolosis/Tb_Test.csv", test, ',')