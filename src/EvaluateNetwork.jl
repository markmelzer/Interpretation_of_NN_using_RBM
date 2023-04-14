using DelimitedFiles
#include("TrainNetwork.jl")
#include("DataPreparation.jl")


params = Hyperparameters(mode = "single")
data = DataParameters(file_name = "./data/HPAIV_train_set_imp_feat.csv", input_length = 112, aa_universe = "ABCDEFGHIKLMNPQRSTVWY?")

#model, avg_acc, avg_loss, AA_dict = train_network_AA(params, data)

m, acc, loss, dict, l, w = train_network_AA(params, data)
val = evaluateTestData(m, "./data/HPAIV_test_set_imp_feat.csv", dict, "./results/output_node_values.csv")

w = m.weight
x = sortperm(abs.(w'[:,1]))
# some kind of cutoff

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

    measures = [0, 0, 0, 0]
    for i in 1:size(encoded_MSA)[1]
        pred = model(encoded_MSA[i,:])
        measures += performance_measure(MSA[i, end], evaluate(pred[1]))
    end
    
    measures
end

# save model
#using BSON: @save
#@save "./results/Exp1_HPAIV/vanilla_model.bson" m


#using BSON: @load
#@load "./results/Exp1_HPAIV/vanilla_model.bson" m