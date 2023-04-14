using Plots
using StatsPlots

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
    ylims!(mini - 0.1, maxi + 0.1)
    #rosetta = [84, 91, 47, 49, 107, 109, 180, 182, 96, 98]
    #vline!(rosetta)
end

function histogram_weights(weights, ylabel="Quantity", xlabel="Weights", bins = 50)
    #=
    Function to plot a histogram of the weights trained in an NN
    =#
    histogram(weights, bins=bins)
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

function heatmap_comparison(rosetta, julia, structure)
    #=
    Funtion to plot a heatmap to compare the found significant features in
    R.ROSETTA and the trained NN (Perceptron)
    =#
    heatmap(hcat(rosetta, julia, structure)', color=palette(:Greys_4), legend=:none)
    xlabel!("Position in AA Sequence")
    ylabel!("Rosetta          Perceptron        Structure")
    title!("Comparison of found significant Features")
end

function plot_losses(losses)
    #=
    Function to plot the loss per epoch.
    =#
    mini = minimum(losses[1])
    maxi = maximum(losses[1])
    x = range(1, length(losses[1]), length(losses[1]))

    plot(x, losses[1])
    if length(losses) > 1
        for loss in losses[2:end]
            plot!(x, loss)
            # mini = minimum([minimum(loss), mini])
            # maxi = maximum([maximum(loss), maxi])
        end
    end

    xlabel!("Epoch")
    xlims!(0, length(loss[1]))
    ylabel!("Loss")
    #ylims!(mini - 0.1, maxi + 0.1)
    #rosetta = [84, 91, 47, 49, 107, 109, 180, 182, 96, 98]
    #vline!(rosetta)
end