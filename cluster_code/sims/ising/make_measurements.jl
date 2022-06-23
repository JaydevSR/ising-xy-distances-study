include("../../src/spinmc.jl")

basepath = "D:/Projects/Dr. Heyl Group/data/"

Temps = [1.8, 1.9, 2.0, 2.1, 2.12, 2.14, 2.16, 2.18, 2.2, 2.22, 2.24, 2.26, 2.28, 2.3, 2.32, 2.34, 2.36, 2.38, 2.4, 2.5, 2.6, 2.7]
τvals = [19 , 19 , 20 , 21 , 21  , 21  , 21  , 21  , 21 , 21  , 21  , 21  , 21  , 21 , 21  , 21  , 21  , 21  , 25 , 30 , 35 , 40 ]

lattice_sizes = [16, 32, 48, 64]
eqsteps = 10000
nconf = 1500

@time ising_getconfigdata_to_txt(lattice_sizes, Temps, nconf, eqsteps, τvals; store_at=basepath, ntau=5)