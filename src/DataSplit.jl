# split data into two training and one test set
using DelimitedFiles
using Random

# set seed to reproduce split
Random.seed!(0)

# read file and check dimensions
file = readdlm("./data/NS1\\NS1.csv", ',')
size(file)

# remove header and id --> unnecessary for training
file = file[2:end,2:end]

# randomly reorder objects 
file= file[shuffle(1:end), :]

# check class distribution for imbalance
sum(file[2:end, end])
size(file)[1] - sum(file[2:end, end]) - 1 # -1 for header 

# split roughly 1500 - 500 - 346: NN probably needs more training data than Rosetta

# get positive and negative objects separately
pos = file[findall(isone, file[2:end, end]).+1,:]
neg = file[findall(iszero, file[2:end, end]).+1,:]

# split positive and negative objects

pos_a = pos[1:347,:]
pos_b = pos[348:462, :]
pos_c = pos[463:end, :]

neg_a = neg[1:1153,:]
neg_b = neg[1154:1537, :]
neg_c = neg[1538:end, :]

# combine positive and negative objects for the three data sets
train1 = vcat(pos_a, neg_a)[shuffle(1:end), :]
train2 = vcat(pos_b, neg_b)[shuffle(1:end), :]
test = vcat(pos_c, neg_c)[shuffle(1:end), :]

# write to files
writedlm("./data/NS1/HPAIV_Train1.csv", train1, ',')
writedlm("./data/NS1/HPAIV_Train2.csv", train2, ',')
writedlm("./data/NS1/HPAIV_Test.csv", test, ',')