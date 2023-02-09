using DelimitedFiles
include("TrainNetwork.jl")
include("DataPreparation.jl")


params = Hyperparameters()
data = DataParameters(file_name = "./data/HPAIV_train_set.csv", input_length = 249, aa_universe = "ABCDEFGHIKLMNPQRSTVWY?")

#model, avg_acc, avg_loss, AA_dict = train_network(params, data)

val = evaluateTestData(model, "./data/HPAIV_test_set.csv", AA_dict, "./results/output_node_values.csv")

function evaluateLayer(model, layer, input_data)
    # assume encoded input data, without the decision (only containing features)
    value = (input_data' * model.layers[layer].weight')' .+ model.layers[layer].bias
    model.layers[layer].σ(value)
end


function evaluateNetwork(model, input_data)
    # return node values for one input
    # assume encoded input data, without the decision (only containing features)
    model_values = [input_data]
    for lay in 1:length(model) 
        #println(model_values, '\t', evaluateLayer(model, lay, input_data))
        append!(model_values, [evaluateLayer(model, lay, input_data)])
    end
    model_values
end


# function to encode data, and run evaluateNetwork for all data points
function evaluateTestData(model, data_file, AA_dict, output_name)
    MSA = readdlm(data_file, ',')
    encoded_MSA = MSA_encode(MSA, AA_dict)
    model_values = [["Output for model with the following amount of layers: ", length(model)]]
    for i in 1:size(encoded_MSA)[1]
        append!(model_values, evaluateNetwork(model, encoded_MSA[i,:]))
    end
    writedlm(output_name, model_values, ',')
    model_values
end