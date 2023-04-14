using Parameters

# create structures to save parameters
@with_kw mutable struct Hyperparameters
    actFunction::Function = Flux.identity
    lossFunction = Flux.logitbinarycrossentropy # attach weights to adjust for class imbalance
    η = 0.0001
    optimizer = Flux.Adam
    epochs = 50
    cv = 5
    seed = 0
    mode = "chain" # single or chain
end

@with_kw mutable struct DataParameters
    file_name::String
    input_length::Int64
    aa_universe::String

end