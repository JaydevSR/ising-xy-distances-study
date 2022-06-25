include("../../src/spinmc.jl")

basepath = "D:/Projects/Dr. Heyl Group/data/"

Temps = readdlm(joinpath(basepath, "ising",  "temperature_values.txt"), ',', Float64)[:, 1]
lattice_sizes = readdlm(joinpath(basepath, "ising", "lattice_sizes.txt"), ',', Int64)[:, 1]

eqsteps = 10000
nconf = 1500
store_mode="w"

for N in lattice_sizes
    @time ising_getconfigdata_to_txt([N], Temps, nconf, eqsteps, τvals; store_at=basepath, ntau=5)
end