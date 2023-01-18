include("SingleLayerPerceptonHPAIV.jl")
include("Statistics.jl")
include("Plots.jl")
using Flux
using CSV

function analysis_nn()
    actFunctions = [Flux.celu, Flux.elu, Flux.gelu, Flux.hardsigmoid, Flux.hardtanh, Flux.leakyrelu, Flux.lisht, Flux.logcosh, Flux.logsigmoid, Flux.mish, Flux.relu, Flux.relu6, Flux.selu, Flux.sigmoid, Flux.softsign, Flux.swish, Flux.identity]
    lossFunctions= [Flux.mse, Flux.msle, Flux.huber_loss, Flux.label_smoothing, Flux.crossentropy, Flux.logitcrossentropy, Flux.binarycrossentropy, Flux.logitbinarycrossentropy, Flux.tversky_loss, Flux.focal_loss]
    ηs = [0.1, 0.01, 0.001, 0.00001]
    optimizers = [Flux.Descent, Flux.Momentum, Flux.Nesterov, Flux.Adam, Flux.AdaMax, Flux.NAdam]
    epochs = [10, 50]
    cvs = [3, 5, 10]
    seeds = [1, 10, 42]
    modes = ["single", "chain"]
    title = ["actFunction", "lossFunction", "η", "optimizer", "epochs", "cv", "seed", "mode", "avg_acc", "avg_loss", "Outliers up", "Outlier down"]
    l = length(actFunctions) * length(lossFunctions) * 4 * 6 * 2 * 9 * 2
    data_tab = [[] for i in 1:l+1]
    data_tab[1] = title
    idx = 1
    for actFunction in actFunctions
        for lossFunction in lossFunctions
            for η in ηs
                for optimizer in optimizers
                    for epoch in epochs
                        for cv in cvs
                            for seed in seeds
                                for mode in modes
                                    params = Hyperparameters(actFunction = actFunction, lossFunction = lossFunction, η = η, optimizer = optimizer, epochs = epoch, cv = cv, seed = seed, mode = mode)
                                    data = DataParameters(file_name = "NS1.csv", input_length = 249, aa_universe = "ABCDEFGHIKLMNPQRSTVWY?")
                                    m, weights, avg_acc, avg_loss = train_perceptron(params, data)
                                    if idx % 1 == 0
                                        println(idx, " of ", l)
                                    end
                                    idx += 1
                                    data_tab[idx] = [actFunction, lossFunction, η, optimizer, epoch, cv, seed, mode, avg_acc, avg_loss, sum(find_outliers(weights)), sum(find_outliers(weights; which = "lower"))]
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    writedlm("analysis.csv", data_tab, ',')
end


