include("../../src/spinmc.jl")

using ScikitLearn

@sk_import neighbors: NearestNeighbors

datapath = "results/xy/xy_config_data.jld2"
println("Using data from: $datapath")
file = jldopen(datapath)

Temps = [0.6, 0.9, 1.2, 1.5]

Nvals = [10, 20, 30, 40]
n_samples = 1000
xcomp=true
ycomp=true

# Nmarks = Dict(10 => :circle, 20 => :rect, 30 => :diamond, 40 => :cross)

for N in Nvals
    f = Figure(resolution = (800, 600))
    ax = Axis(f[1, 1], xlabel = "⟨r₁⟩", ylabel = "Γ", title = "Spin-stiffness vs mean value of r₁ ")

    println("Calculating For N=$(N) ...")
    data = Dict(file["$(N)x$(N)/uncorr_configs"])
    
    for stepT in 1:length(Temps)
        T = Temps[stepT]
        println("   | Temperature = $(T) ...")
        
        if haskey(data, T)
            uncorrelated_spins = data[T]
        else
            println("   |   > No Data Found.")
            continue
        end
    
        uncorrelated_spins = data[T][:, :, 1:n_samples]
    
        println("   |   > Calculating Distances ...")
        # fit the kNN algorithm to uncorrelated spin configs
        model = NearestNeighbors(n_neighbors = 2, algorithm = "ball_tree")
        configs_vec = xy_prepare_vector(uncorrelated_spins; xcomp=xcomp, ycomp=ycomp)
        nnbrs = fit!(model, configs_vec)  # Ignore warning 
    
        # calculate distances
        dists, idxs = NearestNeighbors.kneighbors(nnbrs, configs_vec)
        println("   |   > Done.")
    
        println("   |   > Calculating Spin-stiffness ...")
        configs = [uncorrelated_spins[:, :, i] for i in 1:size(uncorrelated_spins)[3]]
        spin_stiffness = xy_spin_stiffness_config_arr(configs, T, N)
        println("   |   > Done.")

        println("   |   > Plotting Γ v/s r₁ ...")
        scatter!(ax, collect(dists[:, 2]), spin_stiffness, label = "T=$(T)")
        println("   |   > Done.")
    end

    axislegend(ax)

    location = "results/xy/scatter_r1_spin_stiffness_N$(N).pdf"

    println("Saving Plots")
    save(location, f)
    println("Done.")
end

close(file)