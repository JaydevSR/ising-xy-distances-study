using JLD2
using ScikitLearn
using Statistics
using CairoMakie

@sk_import neighbors: NearestNeighbors

lattice_sizes = [10, 20, 30, 40]

# open file for reading configurations
file = jldopen("results/ising/ising_configdata.jld2", "w+")

println("Calculating distances ...")
for N in lattice_sizes
    println("   | For N=$(N)")
    try
        file["$(N)x$(N)/uncorr_configs"]
    catch
        println("   | No Data Found.")
        continue
    end
    uncorr_data = file["$(N)x$(N)/uncorr_configs"]
    
    mean_r1_with_T = []

    for (T, uncorr_spins) in uncorr_data

        # fit the kNN algorithm to uncorrelated spin configs
        model = NearestNeighbors(n_neighbors=2, algorithm="ball_tree")
        nnbrs = fit!(model, [reshape(uncorr_spins[:,:, i], N*N) for i=1:size(uncorr_spins)[3]])  # Ignore warning 
        
        # calculate distances
        dists, idxs = NearestNeighbors.kneighbors(nnbrs, [reshape(uncorr_spins[:,:, i], N*N) for i=1:size(uncorr_spins)[3]])
        push!(mean_r1_with_T, (T, mean(dists[:, 2])))
    end

    # Write generated data to the file
    println("   | Writing data...")
    file["$(N)x$(N)/mean_r1"] = mean_r1_with_T
    println("   | Done.")
end
println("Done.\n")

println("Generting Plots ...")
Temps = first.(file["10x10/mean_r1"])
mean_r1_10 = last.(file["10x10/mean_r1"])
mean_r1_20 = last.(file["20x20/mean_r1"])
f = Figure(resolution = (800, 400), font_size = 12)

ax1 = Axis(f[1,1], xlabel="Temperature", ylabel="⟨r₁⟩", title="Mean value of first NN distance with Temperature")

scatterlines!(ax1, Temps, mean_r1_10, label="Size=10x10", markercolor=:blue)
scatterlines!(ax1, Temps, mean_r1_20, label="Size=20x20", markercolor=:purple)
axislegend()

# display(f)
println("Saving Plot: results/ising/mean_r1_with_T.png")
save("results/ising/mean_r1_with_T.png", f)
println("Done.")

close(file)