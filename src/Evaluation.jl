#=
    This file contains function to evaluate the Neural Network model.
=#

using Parameters

# compute TP, FP, TN, FN
function performance_measure(truth, prediction, pos_value = 1)
    # convert prediction value to Int if returned as array of length 1
    if typeof(prediction) == Vector{Int}
        if length(prediction) == 1
            prediction = prediction[1]
        else
            error("Prediction has length > 1.")
        end
    end
    if truth == prediction
        if truth == pos_value
            return [1, 0, 0, 0]
        else
            return [0, 0, 1, 0]
        end
    else
        if prediction == pos_value
            return [0, 1, 0, 0]
        else
            return [0, 0, 0, 1]
        end
    end
end

# compute loss and accuracy
function loss_and_accuracy(data, model, device, params)
    @unpack lossFunction = params
    acc = 0
    ls = 0.0f0
    num = 0
    measures = [0, 0, 0, 0]
    for (x, y) in data
        x, y = device(x), device(y)
        y_hat = onecold(model(x), 0:1)
        ls += lossFunction(y_hat, y)
        acc += sum(y_hat .== y)
        measures += performance_measure(y, y_hat)
        num += 1
    end
    return ls/num, acc/num, measures
end

# compute sensitivity and specificity
function sensitivity_and_specificity(measures)
    sn = measures[1]/(measures[1] + measures[4])
    sp = measures[3]/(measures[3] + measures[2])
    return sn, sp
end