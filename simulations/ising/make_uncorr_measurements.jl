include("../../src/spinmc.jl")

root_location = "D:/Projects/Dr. Heyl Group/data/ising/"

Temps = [1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0]
τvals = [3  , 4  , 4  , 5  , 5  , 5  , 5  , 5  , 5  , 6  , 7  , 8  , 8  , 8  ]

lattice_sizes = [24, 32, 48, 64]
eqsteps = 10000
n_uncorr = 5000

for N in lattice_sizes
    @time ising_getconfigdata_to_txt([N], Temps, n_uncorr, eqsteps; store_at=root_location, wolff=true)
end