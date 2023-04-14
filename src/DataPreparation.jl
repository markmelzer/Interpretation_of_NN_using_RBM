#=
    This file contains function to prepare data for the use in Neural Networks.
    This includes undersampling, split into test and validation set and encoding of non-numeric data.
=#

using Parameters
using Random

# creating dictionary AA to integers
function AA_to_int(data)
    @unpack aa_universe = data
    nums = collect(1:1:22)./22 #.- 0.5 .- 11 
    shuffle!(nums)
    AA_dict = Dict()
    for i in 1:22
        merge!(AA_dict, Dict(string(aa_universe[i])=>nums[i]))
    end
    AA_dict
end

# function to convert AA in numbers
function encode_AA(AA_dict, aa)
    AA_dict[aa]
end

# function that returns sequence in numbers
function encode_seq(AA_dict, seq)
    encode = zeros(length(seq))
    for i in 1:length(seq)
        encode[i] = AA_dict[seq[i]]
    end
    encode
end

# function that returns MSA in numbers
function MSA_encode(MSA, AA_dict)
    data_encoded = zeros(size(MSA)[1], size(MSA)[2]-1)
    
    for i in 1:size(MSA)[1]
        data_encoded[i,:] = encode_seq(AA_dict, MSA[i, 1:end-1])
    end
    data_encoded
end

# function for undersampling

function undersample(data)
    # get number of objects belonging to smaller decision class
    count = sum(data[:,end])
    len = length(data[:,end])
    println(count, '\t', len)
    println(data[:,end])
    if count > len-count 
        smaller = 0
        bigger = 1
        num = len - count
    else
        smaller = 1
        bigger = 0
        num = count
    end
    # choose same number of bigger decision class randomly
    tmp = data[findall(data[:,end].==bigger),:]
    pos = rand(1:size(tmp)[1], Int(num))

    # return balanced data matrix
    vcat(tmp[pos,:], data[findall(data[:,end].==smaller),:])
end


# get data for amino acid MSA
function get_data_AA(MSA, AA_dict, cross, cv=10)
    # MSA object should not have a header or rownames !
    data_encoded = MSA_encode(MSA, AA_dict)
    # data_encoded = undersample(MSA)
    decision = MSA[:, end]


    upper_test = min(trunc(Int, cross * round(length(decision)/cv)), size(data_encoded)[1])
    lower_test = (cross - 1) * trunc(Int, round(length(decision)/cv)) + 1

    test_data = [(data_encoded[i,:], decision[i]) for i in lower_test:upper_test-1]

    if lower_test == 1
        train_data = [(data_encoded[i,:], decision[i]) for i in upper_test:length(decision)]
    elseif upper_test == length(decision)
        train_data = [(data_encoded[i,:], decision[i]) for i in 1:lower_test-1]
        test_data = vcat(test_data, [(data_encoded[end,:], decision[end])])
    else
        train_data = [(data_encoded[i,:], decision[i]) for i in 1:lower_test-1]
        train_data = vcat(train_data, [(data_encoded[i,:], decision[i]) for i in upper_test:length(decision)])
    end
    
    return train_data, test_data

end


# get data for continuous data
function get_data(MSA)
    # data object should not have a header or rownames, but the decision !
    
    features = MSA[:,1:end-1]
    decision = MSA[1:end, end]

    return [(features[i,:], decision[i]) for i in 1:length(decision)]

end