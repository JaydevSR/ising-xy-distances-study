include("../../src/spinmc.jl")

# open/overwrite file
file = jldopen("results/ising/ising_configdata.jld2", "w+")

# Simulation Parameters
lattice_sizes = [10, 20]
Temps = collect(1.4:0.2:3.6)
corr_times = zeros(Float64, length(Temps))

eqsteps = 2000  # Steps for equilibration
n_uncorr = 50

ising_getconfigdata!(file, lattice_sizes, Temps, n_uncorr, eqsteps; from_infinity=true)

close(file)