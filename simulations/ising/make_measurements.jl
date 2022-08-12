# Make uncorrelated measurements of the Ising model
include("../../src/spinmc.jl")

lattice_sizes = [32, 40, 48, 56, 64]

temperatures = [1.80, 1.90, 2.00, 2.10]
append!(temperatures, [2.20, 2.21, 2.22, 2.23, 2.24, 2.25, 2.26])
append!(temperatures, [2.270, 2.274, 2.278, 2.282, 2.286, 2.290, 2.294, 2.298])
append!(temperatures, [2.302, 2.306, 2.310, 2.314, 2.318, 2.322, 2.33, 2.34, 2.35, 2.36, 2.37, 2.38, 2.39])
append!(temperatures, [2.40, 2.50, 2.60, 2.70])

autocorr_times = [2, 3, 3, 3]
append!(autocorr_times, [6, 6, 7, 9, 10, 13, 16])
append!(autocorr_times, [16, 15, 15, 14, 14, 13, 13, 13])
append!(autocorr_times, [13, 13, 13, 12, 12, 11, 10, 10, 9, 9, 9, 9, 9])
append!(autocorr_times, [8, 12, 15, 24])

nconfigs = 20000
eqsteps = 10000
ntau = 10
store_at = "data/"
start = :cold
verbose = true
mode = "w"
get_mags = true
get_structure_factors = true

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