##
using ScikitLearn

include("../../src/spinmc.jl")

@sk_import neighbors:NearestNeighbors

datapath = "results/xy/xy_config_data.jld2"
println("Using data from: $datapath")
file = jldopen(datapath)

phistr = "phi"
dotstr = "dot"
xstr = "xcomp"
ystr = "ycomp"
emptystr = ""

N = 30
Temps = [0.7, 1.0, 1.3, 1.6]
simtype = "multitemp"
phi = false
xcomp = true
ycomp = true

plotxlab = "r₁ (xcomp=$(xcomp), ycomp=$(ycomp))"
plotylab = "R$(phi ? phistr : dotstr)"

f = Figure(resolution = (800, 600))
ax = Axis(f[1, 1], xlabel = plotxlab, ylabel = plotylab, title = "Lattice size = $(N)")

println("Calculating For N=$(N) ...")

data = Dict(file["$(N)x$(N)/uncorr_configs"])

##

for stepT in 1:length(Temps)
    T = Temps[stepT]
    println("   | Temperature = $(T) ...")
    
    if haskey(data, T)
        uncorrelated_spins = data[T]
    else
        println("   |   > No Data Found.")
        continue
    end

    uncorrelated_spins = data[T]

    println("   |   > Calculating Distances ...")
    # fit the kNN algorithm to uncorrelated spin configs
    model = NearestNeighbors(n_neighbors = 2, algorithm = "ball_tree")
    configs_vec = xy_prepare_vector(uncorrelated_spins; xcomp=xcomp, ycomp=ycomp)
    nnbrs = fit!(model, configs_vec)  # Ignore warning 

    # calculate distances
    dists, idxs = NearestNeighbors.kneighbors(nnbrs, configs_vec)
    println("   |   > Done.")

    mult(a,b) = a*b
    metric = phi ? mult : xy_spindot
    println("   |   > Calculating Structure Factors ...")
    struc_factors = [structure_factor(uncorrelated_spins[:, :, i], N; metric=metric) for i = 1:size(uncorrelated_spins)[3]]
    println("   |   > Done.")

    println("   |   > Plotting R v/s r₁ ...")
    scatter!(ax, collect(dists[:, 2]), struc_factors, label = "T=$(T)")
    println("   |   > Done.")
end

##

r1vals = collect(0:0.01:N*sqrt(2))
Rvals = zeros(Float64, length(r1vals))
@. Rvals = 1 - r1vals^2 / (2*N*N)
lines!(ax, r1vals, Rvals, label="r₁=√(2Nₛ(1-R))", color=:black, linestyle=:dash)

axislegend(ax)
# display(f)
location = "results/xy/$(N)x$(N)/r1$(!ycomp ? xstr : emptystr)$(!xcomp ? ystr : emptystr)_vs_R$(phi ? phistr : dotstr)_$(N)_$(simtype)_xywolff.png"
println("Saving Plots to: $(location)")
save(location, f)
println("Done.")

close(file)
##