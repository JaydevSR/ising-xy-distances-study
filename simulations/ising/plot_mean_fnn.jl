using CairoMakie
using Statistics
using DelimitedFiles
using StaticArrays

basepath = "D:/Projects/Dr. Heyl Group/data/"

Temps = [
    2.00, 2.10, 
    2.11, 2.12, 2.13, 2.14, 2.15, 2.16, 2.17, 2.18, 2.19, 2.20,
    2.21, 2.22, 2.23, 2.24, 2.25, 2.26, 2.27, 2.28, 2.29, 2.30,
    2.31, 2.32, 2.33, 2.34, 2.35, 2.36, 2.37, 2.38, 2.39, 2.40,
    2.50, 2.60
    ]
lattice_sizes = [16, 24, 32, 40, 48, 64]

f = Figure()
ax = Axis(f[1,1], xlabel="temp", ylabel="mean fnn distance")

for N in lattice_sizes
    mean_fnn = zero(Temps)
    Threads.@threads for stepT in eachindex(Temps)
        T = Temps[stepT]
        datafile = joinpath(basepath, "ising", "fnn_dists", "Size$N", "fnn_dists_temp$(T)_size$(N).txt")
        if !isfile(datafile)
            println("| No Data Found (N=$(N), T=$(T)).")
            continue
        end
        fnn_dists = readdlm(datafile, ',', Float64)
        mean_fnn[stepT] = mean(fnn_dists[1:end])
    end
    scatterlines!(Temps, mean_fnn)
end

display(f)