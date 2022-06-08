include("../../src/spinmc.jl")

using ScikitLearn

@sk_import neighbors: NearestNeighbors

basepath = "D:/Projects/Dr. Heyl Group/data/ising/"
println("Using data from: $basepath")

Temps = [2.0, 2.1, 2.2, 2.24, 2.26, 2.265, 2.27, 2.275, 2.28, 2.3, 2.4]

Nvals = [16, 24, 32, 48, 64]
n_samples = 5000

f = Figure(resolution = (1200, 800))
ax = Axis(f[1, 1], xlabel = "T", ylabel = "C (scaled by nsites)", title = "Correlation between r₁ and Ns*R with temperature")
ax2 = Axis(f[1, 2], xlabel = "T", ylabel = "C", title = "Correlation between r₁ and R with temperature")

Nmarks = Dict(16 => :circle, 24 => :rect, 32 => :diamond, 48 => :cross, 64 => :xcross)
Ncols = Dict(16 => :red, 24 => :green, 32 => :blue, 48 => :magenta, 64 => :black)

corr_data = zeros(Float64, (length(Nvals), length(Temps)))

for stepN in 1:length(Nvals)
    N = Nvals[stepN]
    println("Calculating For N=$(N) ...")
    corr_arr = zeros(Float64, length(Temps))
    
    struc_facs_store = Dict([(struc_facs_store[i, 1], struc_facs_store[i, 2:end]) for i=1:size(struc_facs_store)[1]])
    
    Threads.@threads for stepT in eachindex(Temps)
        T = Temps[stepT]
        println("| Process strarted on thread #$(Threads.threadid()): T = $(T)")
        datafile = basepath*"Size$(N)/ising_uncorr_configs_Temp$(Temps[stepT])_N$(N).txt"
        if !isfile(datafile)
            println("| No Measurement Data Found (T=$(T)).")
            println("| Process ended on thread #$(Threads.threadid()): T = $(T)")
            continue
        end
    
        uncorrelated_spins = reshape(readdlm(datafile, ',', Float64), N, N, :)
    
        # fit the kNN algorithm to uncorrelated spin configs
        model = NearestNeighbors(n_neighbors = 2, algorithm = "ball_tree")
        configs_vec = [reshape(uncorrelated_spins[:, :, i], N * N) for i = 1:size(uncorrelated_spins)[3]]
        nnbrs = fit!(model, configs_vec)  # Ignore warning 
    
        # calculate distances
        dists, idxs = NearestNeighbors.kneighbors(nnbrs, configs_vec)
    
        struc_facs_store = basepath*"struc_facs/Size$N/struc_facs_Temps$T_N$(N).txt"
        if !isfile(struc_facs_store)
            println("| No Structure Factor Data Found (T=$(T)).")
            println("| Process ended on thread #$(Threads.threadid()): T = $(T)")
            continue
        end

        struc_facs = readdlm(struc_facs_store, ',', Float64)
        corr_arr[stepT] = mean(dists.*struc_facs) - mean(dists)*mean(struc_facs)
        println("| Process completed on thread #$(Threads.threadid()): T = $(T)")
    end

    println("|")
    println("| Adding Plot ...")
    scatter_points = Point2f.(Temps, corr_arr)
    scatterlines!(ax, scatter_points, color=Ncols[N], label="L=$(N)", marker=Nmarks[N])

    scatter_points2 = Point2f.(Temps, corr_arr/(N^2))
    scatterlines!(ax2, scatter_points2, color=Ncols[N], label="L=$(N)", marker=Nmarks[N])
    println("| Done.")
    println("|")

    corr_data[stepN, 1:length(corr_arr)] = corr_arr
end

axislegend(ax, position = :rb)
axislegend(ax2, position = :rb)

location_plot = "results/ising/r1_R_correlation_plot_v1.pdf"

open(basepath*"r1_R_corr_data.txt", "w") do io
    writedlm(io, corr_data, ',')
end;

open(basepath*"r1_R_corr_Nvals.txt", "w") do io
    writedlm(io, Nvals, ',')
end;

open(basepath*"r1_R_corr_Temps.txt", "w") do io
    writedlm(io, Temps, ',')
end;

println("Saving Plots")
save(location_plot, f)
println("Done.")