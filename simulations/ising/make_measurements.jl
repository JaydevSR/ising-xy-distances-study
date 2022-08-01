# Make uncorrelated measurements of the Ising model
include("../../src/spinmc.jl")

lattice_sizes = [16]#, 32, 48, 64]
# temperatures = [1.80, 1.90, 2.00, 2.10]
# append!(temperatures, [2.20, 2.21, 2.22, 2.23, 2.24, 2.25, 2.26, 2.27, 2.28, 2.29])
# append!(temperatures, [2.30, 2.31, 2.32, 2.33, 2.34, 2.35, 2.36, 2.37, 2.38, 2.39])
# append!(temperatures, [2.40, 2.50, 2.60, 2.70])
temperatures = [2.10]
autocorr_times = [3]
# autocorr_times = [1, 2, 2, 2]
# append!(autocorr_times, [5, 5, 6, 8, 9, 11, 14, 12, 11, 10])
# append!(autocorr_times, [9, 8, 8, 8, 8, 8, 8, 8, 8, 8])
# append!(autocorr_times, [8, 12, 15, 24])

nconfigs = 10000
eqsteps = 10000
ntau = 10
store_at = "data/"
start=:cold
verbose=true
mode="w"
get_mags=true
get_structure_factors=true

@time ising_get_measurements_to_txt(
    lattice_sizes,
    temperatures,
    nconfigs,
    eqsteps,
    autocorr_times;
    ntau=ntau,
    store_at=store_at,
    start=start,
    verbose=verbose,
    mode=mode,
    get_mags=get_mags,
    get_structure_factors=get_structure_factors
    )