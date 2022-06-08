include("../../src/spinmc.jl")

using ScikitLearn

@sk_import neighbors: NearestNeighbors

basepath = "D:/Projects/DQCM (Dr. Heyl)/cluster_data/data/xy/"
println("Using data from: $basepath")

Temps = [0.6, 0.9, 1.2, 1.5]

Nvals = [10, 20, 30, 40]
n_samples = 1000
xcomp=true
ycomp=true

# Nmarks = Dict(10 => :circle, 20 => :rect, 30 => :diamond, 40 => :cross)

for N in Nvals
    f = Figure(resolution = (800, 600))
    ax = Axis(f[1, 1], xlabel = "r₁", ylabel = "Γ", title = "Spin-stiffness vs r₁ (N=$(N))")

    println("Calculating For N=$(N) ...")
    
    for stepT in 1:length(Temps)
        T = Temps[stepT]
        println("   | Temperature = $(T) ...")
        
        datafile = basepath*"Size$(N)/uncorr_configs_Temp$(Temps[stepT])_N$(N).txt"
        
        if !isfile(datafile)
            println("   |   > No Data Found.")
            continue
        end
    
        uncorrelated_spins = reshape(readdlm(datafile, ',', Float64), N, N, :)
    
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
        spin_stiffness = xy_spin_stiffness_config_arr_two(configs, T, N)
        println("   |   > Done.")

        println("   |   > Plotting Γ v/s r₁ ...")
        scatter!(ax, collect(dists[:, 2]), spin_stiffness, label = "T=$(T)")
        println("   |   > Done.")
    end

    axislegend(ax, position=:lt)

    location = "results/xy/scatter_r1_spin_stiffness_term2_N$(N).pdf"

    println("Saving Plots")
    save(location, f)
    println("Done.")
end