using ScikitLearn
using JLD2
using Plots
using Statistics

@sk_import neighbors: NearestNeighbors

lattice_sizes = [10, 20, 30, 40, 50]

# open file for reading configurations
data = jldopen("data/isingdata.jld2", "r")

distance_data = Dict()

println("Calculating distances ...")
for N in lattice_sizes
    println("   | For N=$(N)")
    uncorr_data = data["$(N)x$(N)/uncorr_measurements"]
    nn_dist_with_T = []

    for (T, uncorr_spins) in uncorr_data

        # fit the kNN algorithm to uncorrelated spin configs
        nnbrs = fit!(
            NearestNeighbors(n_neighbors=3, algorithm="ball_tree"),
            # convert to form accepted by algorithm
            [reshape(uncorr_spins[:,:, i], N*N) for i=1:size(uncorr_spins)[3]]
        )  # Ignore warning 
        
        # calculate distances
        dists, idxs = NearestNeighbors.kneighbors(nnbrs, [reshape(uncorr_spins[:,:, i], N*N) for i=1:size(uncorr_spins)[3]])
        push!(nn_dist_with_T, (T, mean(dists[:, 2]), mean(dists[:, 3])))
    end

    distance_data["$(N)"] = nn_dist_with_T
    println("   | Done.")
end
println("Done.\n")

# Plotting
println("Generating Plots ...")
# For 1st NN distance
let N=lattice_sizes[1]
    plot_data = distance_data["$(N)"]
    plot(first.(plot_data), [plot_data[i][2] for i=1:length(plot_data)], label="N=$(N)", markershape=:diamond)
end
for N in lattice_sizes[2:end]
    plot_data = distance_data["$(N)"]
    plot!(first.(plot_data), [plot_data[i][2] for i=1:length(plot_data)], label="N=$(N)", markershape=:diamond)
end
xlabel!("Temperature")
ylabel!("r1")
title!("First NN distance with temperature")
savefig("plots/first_nn_ising.png")

# For 2nd NN distance
let N=lattice_sizes[1]
    plot_data = distance_data["$(N)"]
    plot(first.(plot_data), [plot_data[i][3] for i=1:length(plot_data)], label="N=$(N)", markershape=:diamond)
end
for N in lattice_sizes[2:end]
    plot_data = distance_data["$(N)"]
    plot!(first.(plot_data), [plot_data[i][3] for i=1:length(plot_data)], label="N=$(N)", markershape=:diamond)
end
xlabel!("Temperature")
ylabel!("r2")
title!("Second NN distance with temperature")
savefig("plots/second_nn_ising.png")


# for comparing 1st and 2nd NN distances
let N=10
    plot_data = distance_data["$(N)"]
    plot(first.(plot_data), [plot_data[i][2] for i=1:length(plot_data)], label="r1", markershape=:diamond)
    plot!(first.(plot_data), [plot_data[i][3] for i=1:length(plot_data)], label="r2", markershape=:diamond)
end
xlabel!("Temperature")
ylabel!("distances")
title!("First and Second NN distance with T (N=10)")
savefig("plots/first_second_nn_ising_10.png")

println("Done.")

close(data)