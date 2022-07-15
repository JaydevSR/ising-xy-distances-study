# include("../../src/spinmc.jl")
# using Plots
using CairoMakie

basepath = "D:/Projects/Dr. Heyl Group/data/"

# Temps = readdlm(joinpath(basepath, "ising",  "temperature_values.txt"), ',', Float64)[:, 1]
# lattice_sizes = readdlm(joinpath(basepath, "ising", "lattice_sizes.txt"), ',', Int64)[:, 1]
Temps = [
    2.00, 2.10, 
    2.11, 2.12, 2.13, 2.14, 2.15, 2.16, 2.17, 2.18, 2.19, 2.20,
    2.21, 2.22, 2.23, 2.24, 2.25, 2.26, 2.27, 2.28, 2.29, 2.30,
    2.31, 2.32, 2.33, 2.34, 2.35, 2.36, 2.37, 2.38, 2.39, 2.40,
    2.50, 2.60
    ]
lattice_sizes = [16, 24, 32, 40, 48, 64]

f = Figure(resolution=(1920, 1080))
ax = Axis(f[1, 1], xlabel = "T", ylabel = "C (scaled by nsites)", title = "Correlation between r₁ and Ns*R with temperature")
ax2 = Axis(f[1, 2], xlabel = "T", ylabel = "C", title = "Correlation between r₁ and R with temperature")

Nmarks = Dict(16 => :circle, 24 => :rect, 32 => :diamond, 40 => :rect, 48 => :cross, 64 => :xcross)
Ncols = Dict(16 => :red, 24 => :green, 32 => :blue, 40 => :orange, 48 => :magenta, 64 => :black)

corr_data = zeros(Float64, (length(lattice_sizes), length(Temps)))

for stepN in eachindex(lattice_sizes)
    N = lattice_sizes[stepN]
    println("Calculating For N=$(N) ...")
    corr_arr = zeros(Float64, length(Temps))
    
    for stepT in eachindex(Temps)
        T = Temps[stepT]

        fnn_dists_store = joinpath(basepath, "ising", "fnn_dists", "Size$N", "fnn_dists_temp$(T)_size$(N).txt")
        if !isfile(fnn_dists_store)
            println("| No FNN Distance Data Found at $fnn_dists_store (T=$(T)).")
            println("| Process ended on thread #$(Threads.threadid()): T = $(T)")
            continue
        end
        fnn_dists = readdlm(fnn_dists_store, ',', Float64)[:, 1]
    
        struc_facs_store = joinpath(basepath, "ising", "struc_facs", "Size$N", "struc_facs_temp$(T)_size$(N).txt")
        if !isfile(struc_facs_store)
            println("| No Structure Factor Data Found (T=$(T)).")
            println("| Process ended on thread #$(Threads.threadid()): T = $(T)")
            continue
        end
        struc_facs = readdlm(struc_facs_store, ',', Float64)

        corr_arr[stepT] = mean(fnn_dists.*struc_facs) - mean(fnn_dists)*mean(struc_facs)
        println("| Process completed on thread #$(Threads.threadid()): T = $(T)")
    end

    println("|")
    println("| Adding Plot ...")
    scatterlines!(ax, Temps, corr_arr .* (N^2), color=Ncols[N], label="L=$(N)", marker=Nmarks[N])
    scatterlines!(ax2, Temps, corr_arr, color=Ncols[N], label="L=$(N)", marker=Nmarks[N])
    println("| Done.")
    println("|")

    corr_data[stepN, 1:length(corr_arr)] = corr_arr
end

axislegend(ax, position = :rb)
axislegend(ax2, position = :rb)

location_plot = joinpath("results", "ising", "r1_R_correlation_plot.png")

open(joinpath(basepath, "ising", "r1_R_correlation_data.txt"), "w") do io
    writedlm(io, corr_data, ',')
end;

open(joinpath(basepath, "ising", "temperature_values.txt"), "w") do io
    writedlm(io, Temps, ',')
end;

open(joinpath(basepath, "ising", "lattice_sizes.txt"), "w") do io
    writedlm(io, lattice_sizes, ',')
end;

println("| Saving Plots ...")
save(location_plot, f)
println("| Done.")
println("*===========================================")

display(f)