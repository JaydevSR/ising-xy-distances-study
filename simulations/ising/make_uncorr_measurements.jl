include("../../src/spinmc.jl")

root_location = "D:/Projects/Dr. Heyl Group/data/ising/"

Temps = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]

lattice_sizes = [16, 24, 32, 48, 64]
eqsteps = 10000
n_uncorr = 10000

ising_getconfigdata_to_txt(lattice_sizes, Temps, n_uncorr, eqsteps; store_at=root_location, calculate_autocorr_times=true, wolff=true)