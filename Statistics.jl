using Statistics


function find_whisker(weights; which="upper", α=1.5)
    #=
    Function that finds the value of a whisker of a boxplot
    =#
    IQR = quantile(weights, 0.75) - quantile(weights, 0.25)
    if which == "upper"
        return quantile(weights, 0.75) + α*IQR
    elseif which == "lower"
        return quantile(weights, 0.25) - α*IQR
    else
        error("find_whisker: argument 'which' must be 'upper' or 'lower'")
    end
end

function find_outliers(data; which="upper", α=1.5)
    #=
    Function that finds the outliers of a boxplot.
    Either above or below it.
    =#
    whisker = find_whisker(data, which, α)
    if which == "upper"
        return data .> whisker
    elseif which == "lower"
        return data .< whisker
    else
        error("find_whisker: argument 'which' must be 'upper' or 'lower'")
    end
end

function preprocess_rosetta(rosetta, numFeatures)
    #=
    Function that preprocesses an input array given from ROSETTA.
    Input array: List of lists of significant features
    E.g., [[Rules of length 1], [Rules of length 2], ...]
    Output: Array of length numFeatures, that indicates the length of the rule a
        feature belongs to. In case of multiple rule, the one with the shortest
        length is choose.
    =#
    out = zeros(Int16, numFeatures)
    for i in length(rosetta):-1:1
        for j in rosetta[i]
            out[j] = i
        end
    end
    out
end