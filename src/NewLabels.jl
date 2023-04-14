using DelimitedFiles

file_name = "./data/HPAIV_test_set_imp_feat.csv"

secTrain = readdlm(file_name, ',')

# function to encode data, and run evaluateNetwork for all data points
function evaluateTestData(model, MSA, AA_dict)
    encoded_MSA = MSA_encode(MSA, AA_dict)

    measures = [0, 0, 0, 0]
    p = []
    for i in 1:size(encoded_MSA)[1]
        pred = model(encoded_MSA[i,:])
        append!(p, pred)
        measures += performance_measure(MSA[i, end], evaluate(pred[1]))
    end
    
    return(p, measures)
end

pred, meas = evaluateTestData(m, secTrain, dict)
newLabel = evaluate.(pred)

newData = hcat(secTrain, newLabel)
writedlm("AADataNewLabel.csv", newData, ',')