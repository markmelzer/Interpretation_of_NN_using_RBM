# split data into two training and one test set
using DelimitedFiles

file = readdlm("./data/NS1\\NS1.csv", ',')
size(file)

file = file[:,2:end]
sum(file[2:end, end])
size(file)[1] - sum(file[2:end, end]) - 1 # -1 for header 

# split roughly 1500 - 500 - 346: NN probably needs more training data than Rosetta
pos = file[findall(isone, file[2:end, end]).+1,:]
neg = file[findall(iszero, file[2:end, end]).+1,:]

a = 1500/2346
b = 500/2346
c = 346/2346

a*1803
b*1803
c*1803

pos_a = pos[1:347,:]
pos_b = pos[348:462, :]
pos_c = pos[463:end, :]

neg_a = neg[1:1153,:]
neg_b = neg[1154:1537, :]
neg_c = neg[1538:end, :]

train1 = vcat(pos_a, neg_a)[shuffle(1:end), :]
train2 = vcat(pos_b, neg_b)[shuffle(1:end), :]
test = vcat(pos_c, neg_c)[shuffle(1:end), :]

writedlm("./data/HPAIV_Train1.csv", train1, ',')