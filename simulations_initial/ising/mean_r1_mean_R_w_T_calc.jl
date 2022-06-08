include("../../src/spinmc.jl")

using ScikitLearn

@sk_import neighbors: NearestNeighbors

use_wolff = false
wolff_str = "isingwolff"
metro_str = "isingmetro"

Temps = collect(1.8:0.1:2.6)

N = 20
spins = rand([-1.0, 1.0], (N, N))
eqsteps = 2000
n_uncorr = 100
println("Calculating for N=$(N) ...")

mean_r1 = zeros(Float64, length(Temps))
mean_R = zeros(Float64, length(Temps))

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

    mean_r1[stepT] = mean(dists[:, 2])
    mean_R[stepT] = mean(struc_factors)
end

##
println("Plotting ⟨R⟩ v/s ⟨r₁⟩ ...")
f = Figure(resolution = (800, 600))
ax = Axis(f[1, 1], xlabel = "⟨r₁⟩", ylabel = "⟨R⟩", title = "Lattice size = $(N)")
# cols = cgrad(:matter, length(Temps), categorical=true, rev=true)
scatter_points = Point2f.(mean_r1, mean_R)
plt = scatter!(ax, scatter_points, color=1:length(Temps))
Colorbar(
    f[1,2],
    # colormap = cgrad(:matter, rev=true),
    limits = (Temps[1], Temps[end]),
    ticks = Temps[1:end], tickalign=1,
    label="Temperature"
)
# display(f)
save("results/ising/mean_r1_mean_R_w_T_$(N)_$(use_wolff ? wolff_str : metro_str).png", f)
println("Done.")
##