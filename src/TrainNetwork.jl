using Flux
using Random
using Flux.Data: DataLoader
using Flux: onehotbatch, onecold, onehot, @epochs
using NNlib
using DelimitedFiles
include("DataPreparation.jl")
include("Evaluation.jl")
include("ModelDefinitions.jl")
include("Parameters.jl")


# call parameters
params = Hyperparameters()
data = DataParameters(file_name = "NS1.csv", input_length = 10, aa_universe = "ABCDEFGHIKLMNPQRSTVWY?")




function train_perceptron(params, data)
    device = cpu
    @unpack lossFunction, η, optimizer, epochs, cv, seed, mode = params
    @unpack file_name, input_length = data



    Random.seed!(seed)

    # get data
    AA_dict = AA_to_int(data)
    MSA = readdlm(file_name, ',')
    
    # undersampling might be necessary
    #MSA = undersample(MSA[1:end,:])

    #shuffle MSA (very important, otherwise very bad results)
    MSA = MSA[shuffle(1:end), :]

    #idx = [90, 86, 89, 87, 88, 229, 224, 49, 138, 182, 210, 129, 219, 177, 64, 218, 163, 237, 233, 68, 56, 238, 164, 156, 109, 60, 101, 98, 217, 114, 209, 61, 72, 71, 82, 112, 28, 181, 223, 140, 97, 26, 122, 174, 7, 57, 127, 148, 203, 19, 150, 29, 15, 25, 34, 204, 106, 206, 27, 24, 55, 239, 157, 235, 221, 231, 85, 240, 232, 149, 22, 123, 169, 242, 45, 236, 128, 43, 125, 119, 105, 241, 8, 227, 234, 23, 192, 74, 216, 225, 117, 135, 230, 197, 154, 102, 139, 151, 201, 63, 251]  
    #MSA = MSA[:, vcat(idx[1:input_length], idx[end])]

    #MSA[findall(MSA[:,end].==0),end] .= -1
    #MSA[findall(MSA[:,end].==1),end] .= 1
    

    avg_acc = 0
    avg_loss = 0

    model = 0
    weights = zeros(input_length)

    # TODO: make functions for epoch, optimizing part etc

    for cross in 1:cv
        train_data, test_data = get_data(MSA, AA_dict, cross, cv)

        # construct model
        model = create_perceptron(params, data) |> device

        # optimizer
        opt = Flux.setup(optimizer(η), model)
        
        # Training
        test_loss, test_acc = 0, 0
        for epoch in 1:epochs
            # shuffle!(train_data)
            for (x, y) in train_data
                #x, y = device(x), device(y)
                gs = gradient(m -> lossFunction(m(x), y), model)
                Flux.Optimise.update!(opt, model, gs[1])
            end
        end

        if mode == "single"
            weights .+= abs.(model.weight[1,:]/cv)
            #weights .+= (model.weight[1,:]/cv)
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
    end
    println("Avg. accuracy: ", avg_acc, "\t avg. loss: $avg_loss", "\t pos. predicted: $pos")
    
    # return last model
    return model, weights, avg_acc, avg_loss, hcat(MSA_encode(MSA, AA_dict), MSA[:,end]), MSA
end

