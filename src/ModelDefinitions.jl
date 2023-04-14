using Parameters
using Flux


# function that makes Perceptron
function create_perceptron(params, data)
    @unpack input_length = data
    @unpack actFunction, mode = params
    #Chain(
    #    Flux.Dense(input_length => 10, relu),
    #    Flux.Dense(10=>2, x->σ.(x))
    #)
    if mode == "chain"
        return Flux.Chain(Flux.Dense(input_length => input_length, actFunction),
        Flux.Dense(input_length => input_length, actFunction),
        Flux.Dense(input_length => input_length, actFunction),
        Flux.Dense(input_length => input_length, actFunction),
        Flux.Dense(input_length => input_length, actFunction),
        Flux.Dense(input_length => input_length, actFunction),
        Flux.Dense(input_length => input_length, actFunction),
        Flux.Dense(input_length => input_length, actFunction),
        Flux.Dense(input_length => input_length, actFunction),
        Flux.Dense(input_length => 1, actFunction))
    end
    Flux.Dense(input_length => 1, actFunction)
end