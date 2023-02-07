include("SingleLayerPerceptonHPAIV.jl")

params = Hyperparameters()
data = DataParameters(file_name = "NS1.csv", input_length = 249, aa_universe = "ABCDEFGHIKLMNPQRSTVWY?")

model, avg_weights, avg_acc, avg_loss, encoded_data, data = train_perceptron(params, data)

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