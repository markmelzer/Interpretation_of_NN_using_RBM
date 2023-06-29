using Parameters
using Flux


# function that makes Perceptron
function create_perceptron(params, data)
    @unpack input_length = data
    @unpack actFunction, mode = params
    input_length *= 7
    if mode == "chain"
        return Flux.Chain(Flux.Dense(input_length => input_length, actFunction),
        Flux.Dense(input_length => input_length, actFunction),
        Flux.Dense(input_length => input_length, actFunction),
        Flux.Dense(input_length => 2, actFunction))
    end
    Flux.Chain(Flux.Dense(input_length => 2, actFunction), softmax)
end