include("TrainNetwork.jl")


params = Hyperparameters()
data = DataParameters(file_name = "./data/HPAIV_train_set.csv", input_length = 249, aa_universe = "ABCDEFGHIKLMNPQRSTVWY?")

model, avg_acc, avg_loss, encoded_data, data = train_network(params, data)

l = 1

layer = model.layers[l]
weights = layer.weight
bias = layer.bias

# TODO: this is only for linear activation function in this layer!!! (can check it out with layer.σ)
values_in_layer = (encoded_data[:,1:end-1] * weights')' .+ bias


# write to files: encoded amino acid sequence + decision, values in layers
using DelimitedFiles
writedlm("UsedMSANS1.csv", data, ',')
writedlm("NodeValues2L3N.csv", values_in_layer, ',')


# return model

# for each layer 
    # compute values of nodes, based on input (or values of previous network), weights, bias and activation function