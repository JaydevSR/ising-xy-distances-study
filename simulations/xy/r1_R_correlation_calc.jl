include("../../src/spinmc.jl")

using ScikitLearn

@sk_import neighbors: NearestNeighbors

datapath = "results/xy/xy_config_data.jld2"
println("Using data from: $datapath")
file = jldopen(datapath)

Temps = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]

Nvals = [10, 20, 30, 40]
n_samples = 500
xcomp=true
ycomp=true

f = Figure(resolution = (800, 600))
ax = Axis(f[1, 1], xlabel = "T", ylabel = "⟨R⋅r₁⟩ - ⟨R⟩⟨r₁⟩ ", title = "Correlation between r₁ and R with temperature")

Nmarks = Dict(10 => :circle, 20 => :rect, 30 => :diamond, 40 => :cross)
Ncols = Dict(10 => :red, 20 => :green, 30 => :blue, 40 => :black)

for N in Nvals
    println("Calculating For N=$(N) ...")
    data = Dict(file["$(N)x$(N)/uncorr_configs"])

    corr_arr = zeros(Float64, length(Temps))
    
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
    
        println("   |   > Calculating Structure Factors ...")
        struc_factors = [structure_factor(uncorrelated_spins[:, :, i], N; metric=xy_spindot, scaled=true) for i = 1:size(uncorrelated_spins)[3]]
        println("   |   > Done.")

        corr_arr[stepT] = mean(dists.*struc_factors) - mean(dists)*mean(struc_factors)
    end

    println("   |   > Adding Plot ...")
    scatter_points = Point2f.(Temps, corr_arr)
    scatterlines!(ax, scatter_points, color=Ncols[N], label="L=$(N)", marker=Nmarks[N])
    println("   |   > Done.")
end

axislegend(ax, position = :rb)

location1 = "results/xy/r1_R_correlation.pdf"

println("Saving Plots")
save(location1, f)
println("Done.")

close(file)