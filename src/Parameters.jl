using Parameters

# create structures to save parameters
@with_kw mutable struct Hyperparameters
    actFunction::Function = Flux.celu
    lossFunction = Flux.logitbinarycrossentropy
    η = 0.0001
    optimizer = Flux.Adam
    epochs = 10
    cv = 10
    seed = 42
    mode = "chain" # single or chain
end

@with_kw mutable struct DataParameters
    file_name::String
    input_length::Int64
    aa_universe::String

end