include("../../src/spinmc.jl")

file = jldopen("results/xy/xy_config_data.jld2", "w+")

Temps = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]
τvals = [20 , 25 , 30 , 35 , 45 , 55 , 65 , 75 , 85 , 95 , 115, 145]

lattice_sizes = [40] # [10, 20, 30]
eqsteps = 2000
n_uncorr = 1000

xy_getconfigdata!(file, lattice_sizes, Temps, n_uncorr, eqsteps; autocorr_times = τvals)

close(file)