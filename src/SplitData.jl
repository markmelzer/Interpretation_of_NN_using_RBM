using StatsBase
using Random

# set seed
Random.seed!(6)

# read file, remove Protein name and position
file = readdlm("NS1.csv", ',')
file = file[:,2:end]
file = file[2:end,:]

# get number of sequences
num = size(file)[1]

# split data into pathogenic and non-pathogenic sequences
patho = file[file[:,end] .== 1, :]
non_patho = file[file[:,end] .== 0, :]

# select 10% of data for testing (with equal amount of pathogenic and non-pathogenic sequences)
pos1 = sample(1:size(patho)[1], convert(Int, floor(num/20)), replace = false)
pos2 = sample(1:size(non_patho)[1], convert(Int, floor(num/20)), replace = false)

# concat test and train set, then shuffle
test = vcat(patho[pos1,:], non_patho[pos1,:])
train = vcat(patho[1:end .∉ [pos2],:], non_patho[1:end .∉ [pos2],:])

test = test[shuffle(1:end),:]
train = train[shuffle(1:end),:]

# save into files
writedlm("HPAIV_test_set.csv", test, ',')
writedlm("HPAIV_train_set.csv", train, ',')