using Flux
using Random
using Flux.Data: DataLoader
using Flux: onehotbatch, onecold, onehot, @epochs
using NNlib
using DelimitedFiles
using Parameters



# create structures to save parameters
@with_kw mutable struct Hyperparameters
    actFunction::Function = Flux.identity
    lossFunction = Flux.mse
    η = 0.0001
    optimizer = Flux.Adam
    epochs = 10
    cv = 10
    seed = 10
    mode = "single" # single or chain
end

@with_kw mutable struct DataParameters
    file_name::String
    input_length::Int64
    aa_universe::String

end
# create structures
params = Hyperparameters()
data = DataParameters(file_name = "NS1.csv", input_length = 249, aa_universe = "ABCDEFGHIKLMNPQRSTVWY?")

# creating dictionary AA to integers
function AA_to_int(data)
    @unpack aa_universe = data
    nums = collect(1:1:22)./22 #.- 0.5 .- 11 
    shuffle!(nums)
    AA_dict = Dict()
    for i in 1:22
        merge!(AA_dict, Dict(string(aa_universe[i])=>nums[i]))
    end
    println(AA_dict)
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
    data_encoded = zeros(size(MSA)[1]-1, size(MSA)[2]-2)
    
    for i in 1:size(MSA)[1]-1
        data_encoded[i,:] = encode_seq(AA_dict, MSA[i+1, 2:end-1])
    end
    data_encoded
end

# function for undersampling

function undersample(data)
    # get number of objects belonging to smaller decision class
    count = sum(data[:,end])
    len = length(data[:,end])
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
    pos = rand(1:size(tmp)[1], num)

    # return balanced data matrix
    vcat(tmp[pos,:], data[findall(data[:,end].==smaller),:])
end


# get data 
function get_data(MSA, AA_dict, cross, cv=10)
    data_encoded = MSA_encode(MSA, AA_dict)
    decision = MSA[2:end, end]
    # split data to test and train
    x = trunc(Int, length(decision)/cv)
    lower_test = (cross - 1) * x + 1
    upper_test = cross * x

    test_data = [(data_encoded[i,:], decision[i]) for i in lower_test+1:upper_test]

    if lower_test == 1
        train_data = [(data_encoded[i,:], decision[i]) for i in upper_test+1:length(decision)]
    elseif upper_test == length(decision)
        train_data = [(data_encoded[i,:], decision[i]) for i in 1:lower_test]
    else
        train_data = [(data_encoded[i,:], decision[i]) for i in 1:lower_test]
        train_data = vcat(train_data, [(data_encoded[i,:], decision[i]) for i in upper_test+1:length(decision)])
    end

    return train_data, test_data

end

# function that makes Perceptron
function create_perceptron(params, data)
    @unpack input_length = data
    @unpack actFunction, mode = params
    #Chain(
    #    Flux.Dense(input_length => 10, relu),
    #    Flux.Dense(10=>2, x->σ.(x))
    #)
    if mode == "chain"
        return Flux.Chain(Flux.Dense(input_length => 1, actFunction), make_binary)
    end
    Flux.Dense(input_length => 1, actFunction)
end

# function to make classification binary
function make_binary(values)
    #println(values, evaluate.(values))
    evaluate.(values)
end

# function to assign value to classification
function evaluate(y_hat)
    if y_hat > 0
        return 1
    end
    return -1

end

function loss_and_accuracy(data, model, device, params)
    @unpack lossFunction = params
    acc = 0
    ls = 0.0f0
    num = 0
    pos_pred = 0
    for (x, y) in data
        x, y = device(x), device(y)
        y_hat = model(x)
        #ls += lossFunction(evaluate.(y_hat), y, agg=sum)
        pos_pred += sum(evaluate.(y_hat).==1)
        ls += lossFunction(evaluate.(y_hat)[1], y)
        acc += sum(evaluate.(y_hat) .== y)
        #println(y_hat[1], '\t', evaluate.(y_hat), '\t', y)
        #ls += loss(y_hat[1], y, agg=sum)
        #acc += sum(y_hat[1] .== y)
        num += 1
    end
    return ls/num, acc/num, pos_pred/num
end

function train_perceptron(params, data)
    device = cpu
    @unpack lossFunction, η, optimizer, epochs, cv, seed, mode = params
    @unpack file_name = data



    Random.seed!(seed)

    # get data
    AA_dict = AA_to_int(data)
    MSA = readdlm(file_name, ',')

    #shuffle MSA (very important, otherwise very bad results)
    MSA = undersample(MSA[2:end,:])
    MSA = MSA[shuffle(2:end), :]

    MSA[findall(MSA[:,end].==0),end] .= -1
    

    avg_acc = 0
    avg_loss = 0

    model = 0
    weights = zeros(249)

    pos = 0

    # TODO: make functions for epoch, optimizing part etc

    for cross in 1:cv
        train_data, test_data = get_data(MSA, AA_dict, cross, cv)
        # construct model
        model = create_perceptron(params, data) |> device # 249: length of MSA

        # optimizer
        opt = Flux.setup(optimizer(η), model)
        
        # Training
        test_loss, test_acc = 0, 0
        for epoch in 1:epochs
            # shuffle!(train_data)
            for (x, y) in train_data
                #x, y = device(x), device(y)
                #TODO: trains for 0 and 1, not positive or negative 
                gs = gradient(m -> lossFunction(m(x), y), model)
                Flux.Optimise.update!(opt, model, gs[1])
            end
        end

        if mode == "single"
            weights .+= abs.(model.weight[1,:]/cv)
        elseif mode == "chain"
            weights .+= abs.(model.layers[1].weight[1,:]/cv)
        else
            println("Illegal model choice")
        end
        #train_loss, train_acc = loss_and_accuracy(train_data, model, device)
        test_loss, test_acc, pos_pred = loss_and_accuracy(test_data, model, device, params)
        #println(" train_loss = $train_loss, train_accuracy = $train_acc")
        println(" test_loss = $test_loss, test_accuracy = $test_acc, pos. predicted: $pos_pred")
        
        avg_acc += test_acc/cv
        avg_loss += test_loss/cv
        pos += pos_pred/cv
    end
    println("Avg. accuracy: ", avg_acc, "\t avg. loss: $avg_loss", "\t pos. predicted: $pos")
    
    # return last model
    return model, weights, avg_acc, avg_loss
end

