include("../../src/spinmc.jl")

# file = jldopen("results/xy/xy_config_data.jld2", "w+"; compress=true)

# Temps = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]
# τvals = [20 , 25 , 30 , 35 , 45 , 55 , 65 , 75 , 85 , 95 , 115, 145]

# lattice_sizes = [10, 20, 30, 40]
# eqsteps = 2000
# n_uncorr = 1000

# xy_getconfigdata!(file, lattice_sizes, Temps, n_uncorr, eqsteps; autocorr_times = τvals)

# close(file)

root_location = "results/xy/measurements/"

Temps = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]
τvals = [20 , 25 , 30 , 35 , 45 , 55 , 65 , 75 , 85 , 95 , 115, 145]

lattice_sizes = [10]
eqsteps = 2000
n_uncorr = 1000

xy_getconfigdata_to_txt(lattice_sizes, Temps, n_uncorr, eqsteps; store_at=root_location, autocorr_times = τvals)