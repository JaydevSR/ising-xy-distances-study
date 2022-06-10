using CairoMakie
using Statistics
using DelimitedFiles

Temps = [1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0]

Nvals = [16, 24, 32, 48]

f = Figure()
ax = Axis(f[1,1], xlabel="temp", ylabel="mean fnn distance")

for N in Nvals
    mean_fnn = zero(Temps)
    for stepT in eachindex(Temps)
        T = Temps[stepT]
        datafile = basepath*"fnn_dists/Size$N/fnn_dists_Temp$(T)_N$(N).txt"
        if !isfile(datafile)
            println("| No Data Found (N=$(N), T=$(T)).")
            continue
        end
        fnn_dists = readdlm(datafile, ',', Float64)
        mean_fnn[stepT] = mean(fnn_dists[1:100])
    end
    scatterlines!(Temps, mean_fnn)
end

display(f)