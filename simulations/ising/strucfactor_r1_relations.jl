using ScikitLearn

include("../../src/spinmc.jl")

@sk_import neighbors:NearestNeighbors

N = 30
Temps = [1.8, 2.2, 2.6]

use_wolff = true
wolff_str = "isingwolff"
metro_str = "isingmetro"

f = Figure(resolution = (800, 600))
ax = Axis(f[1, 1], xlabel = "r₁", ylabel = "R", title = "Lattice size = $(N)")

println("Calculating For N=$(N) ...")

spins = rand([-1.0, 1.0], (N, N))
eqsteps = 2000
n_uncorr = 1000

for stepT in 1:length(Temps)
    T = Temps[stepT]
    println("   | Temperature = $(T) ...")

    ising_equilibrate_system!(spins, T, eqsteps)

    println("   |   > Calculating correlation time ...")
    τ = ising_getcorrtime!(spins, T; wolff = use_wolff)
    println("   |   > Done.")

    println("   |   > Making uncorrelated measurements (τ=$(τ)) ...")
    uncorrelated_spins = ising_getuncorrconfigs!(spins, T, τ, n_uncorr; wolff = use_wolff)
    println("   |   > Done.")

    println("   |   > Calculating Distances ...")
    # fit the kNN algorithm to uncorrelated spin configs
    model = NearestNeighbors(n_neighbors = 2, algorithm = "ball_tree")
    nnbrs = fit!(model, [reshape(uncorrelated_spins[:, :, i], N * N) for i = 1:size(uncorrelated_spins)[3]])  # Ignore warning 

    # calculate distances
    dists, idxs = NearestNeighbors.kneighbors(nnbrs, [reshape(uncorrelated_spins[:, :, i], N * N) for i = 1:size(uncorrelated_spins)[3]])
    println("   |   > Done.")

    println("   |   > Calculating Structure Factors ...")
    struc_factors = [structure_factor(uncorrelated_spins[:, :, i], N) for i = 1:size(uncorrelated_spins)[3]]
    println("   |   > Done.")

    println("   |   > Plotting R v/s r₁ ...")
    scatter!(ax, collect(dists[:, 2]), struc_factors, label = "T=$(T)")
    println("   |   > Done.")
end

r1vals = collect(0:0.01:N*sqrt(2))
Rvals = zeros(Float64, length(r1vals))
@. Rvals = 1 - r1vals^2 / (2*N*N)
lines!(ax, r1vals, Rvals, label="r₁=√(2Nₛ(1-R))", color=:black, linestyle=:dash)

axislegend(ax)
display(f)
println("Saving Plots to: results/ising/r1_vs_R_$(N)_$(use_wolff ? wolff_str : ising_str).png")
save("results/ising/r1_vs_R_$(N)_$(use_wolff ? wolff_str : metro_str).png", f)
println("Done.")