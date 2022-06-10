include("../../src/spinmc.jl")

root_location = "D:/Projects/Dr. Heyl Group/data/ising/"

Temps = [1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0]
# τvals = [3  , 4  , 4  , 5  , 5  , 5  , 5  , 5  , 5  , 6  , 7  , 8  , 8  , 8  ]
τvals = [3  , 4  , 8  , 11 , 20 , 185, 200, 121, 42 , 34 , 26 , 20 , 21 , 16 ]

lattice_sizes = [16, 24, 32, 48]
eqsteps = 5000
n_uncorr = 2000

for N in lattice_sizes
    @time ising_getconfigdata_to_txt([N], Temps, n_uncorr, eqsteps, τvals; store_at=root_location, wolff=false, ntau=3)
end