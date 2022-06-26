using CairoMakie
using Statistics
using DelimitedFiles

basepath = "D:/Projects/Dr. Heyl Group/data/"

Temps = readdlm(joinpath(basepath, "ising",  "temperature_values.txt"), ',', Float64)[:, 1]
lattice_sizes = readdlm(joinpath(basepath, "ising", "lattice_sizes.txt"), ',', Int64)[:, 1]

f = Figure()
ax = Axis(f[1,1], xlabel="temp", ylabel="mean fnn distance")

for N in lattice_sizes
    mean_fnn = zero(Temps)
    for stepT in eachindex(Temps)
        T = Temps[stepT]
        # datafile = basepath*"ising/struc_facs/Size$N/struc_facs_temp$(T)_size$(N).txt"
        datafile = basepath*"ising/fnn_dists/Size$N/fnn_dists_temp$(T)_size$(N).txt"
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