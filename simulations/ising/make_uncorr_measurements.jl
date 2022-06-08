include("../../src/spinmc.jl")

root_location = "D:/Projects/Dr. Heyl Group/data/ising/"

Temps = [2.0, 2.1, 2.2, 2.24, 2.26, 2.265, 2.27, 2.275, 2.28, 2.3, 2.4]
τvals = [5  , 5  , 5  , 5   , 5   , 5    , 5   , 5    , 5   , 5  , 5  ]

lattice_sizes = [16, 24, 32, 48, 64]
eqsteps = 10000
n_uncorr = 5000

for N in lattice_sizes
    @time ising_getconfigdata_to_txt([N], Temps, n_uncorr, eqsteps; store_at=root_location, wolff=true)
end