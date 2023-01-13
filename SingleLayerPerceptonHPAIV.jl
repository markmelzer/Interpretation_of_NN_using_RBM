using Flux
using Random
using Flux.Data: DataLoader
using Flux: onehotbatch, onecold, onehot, @epochs

# read data
using DelimitedFiles

# creating dictionary AA to integers
function AA_to_int()
    nums = collect(1:1:21)
    shuffle!(nums)
    AA_dict = Dict()
    aa_universe = "ACDEFGHIKLMNPQRSTVWY?"
    for i in 1:21
        merge!(AA_dict, Dict(aa_universe[i]=>nums[i]))
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
        encode[i] = AA_dict[only(seq[i])]
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
function create_perceptron(input_length)
    #Chain(
    #    Flux.Dense(input_length => 10, relu),
    #    Flux.Dense(10=>2, x->σ.(x))
    #)
    Flux.Dense(input_length => 1, identity)
end

# function to assign value to classification
function evaluate(y_hat)
    if y_hat > 0
        return 1
    end
    return 0
    
end



# function that trains Perceptron on data
const loss = Flux.logitbinarycrossentropy

function loss_and_accuracy(data, model, device)
    acc = 0
    ls = 0.0f0
    num = 0
    for (x, y) in data
        x, y = device(x), device(y)
        y_hat = model(x)
        ls += loss(evaluate.(y_hat), y, agg=sum)
        acc += sum(evaluate.(y_hat) .== y)
        num += 1
    end
    return ls/num, acc/num
end

function train_perceptron(file::String, cv=10, epochs = 7, seed = 42)
    device = cpu
    Random.seed!(seed)

    # get data
    AA_dict = AA_to_int()
    MSA = readdlm(file, '\t')

    avg_acc = 0
    avg_loss = 0

    model = 0
    weights = zeros(249)
    for cross in 1:cv
        train_data, test_data = get_data(MSA, AA_dict, cross, cv)
        # construct model
        model = create_perceptron(249) |> device # 249: length of MSA

        # optimizer
        η::Float64 = 0.1
        opt = Flux.setup(Descent(η), model)
        
        # Training
        test_loss, test_acc = 0, 0
        for epoch in 1:epochs
            # shuffle!(train_data)
            for (x, y) in train_data
                #x, y = device(x), device(y)
                #TODO: trains for 0 and 1, not positive or negative 
                gs = gradient(m -> loss(m(x), y), model)
                Flux.Optimise.update!(opt, model, gs[1])
            end
        end


        weights .+= abs.(model.weight[1,:]/cv)
        # report on train and test
        #train_loss, train_acc = loss_and_accuracy(train_data, model, device)
        test_loss, test_acc = loss_and_accuracy(test_data, model, device)
        #println(" train_loss = $train_loss, train_accuracy = $train_acc")
        println(" test_loss = $test_loss, test_accuracy = $test_acc")
        
        avg_acc += test_acc/cv
        avg_loss += test_loss/cv
    end
    println("Avg. accuracy: ", avg_acc, "\t avg. loss: $avg_loss")
    
    # return last model
    return model, weights
end

# function that evaluates performance

# function that plots weights over position (check for bias)