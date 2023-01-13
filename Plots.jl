using Plots

function plot_weights(weights)
    #=
    Function to plot the weights trained in an NN over the positions
    of an Amino Acid sequence.
    Vertical Lines are added to highlight the regions found most
    significant by R.ROSETTA 
    =#
    mini = minimum(weights)
    maxi = maximum(weights)
    x = range(1, length(weights), length(weights))
    plot(x, weights)
    xlabel!("AA Position")
    xlims!(0, length(weights))
    ylabel!("Weight")
    ylims!(mini - 1, maxi + 1)
    rosetta = [84, 91, 47, 49, 107, 109, 180, 182, 96, 98]
    vline!(rosetta)
end

function histogram_weights(weights, ylabel="Quantity", xlabel="Weights")
    #=
    Function to plot a histogram of the weights trained in an NN
    =#
    histogram(weights)
    ylabel!(ylabel)
    xlabel!(xlabel)
    plot!(legend=false)
end

function boxplot_weights(weights, ylabel="Weight")
    #=
    Function to plot a boxplot of the weights trained in an NN
    =#
    boxplot(weights)
    ylabel!(ylabel)
    plot!(legend=false)
end

function heatmap_comparison(rosetta, julia)
    #=
    Funtion to plot a heatmap to compare the found significant features in
    R.ROSETTA and the trained NN (Perceptron)
    =#
    heatmap(hcat(rosetta, julia)', color=cgrad(:binary), legend=:none)
    xlabel!("Position in AA Sequence")
    ylabel!("Rosetta                  Perceptron")
    title!("Comparison of found significant Features")
end