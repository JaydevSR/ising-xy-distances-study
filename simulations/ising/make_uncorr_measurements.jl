include("../../src/spinmc.jl")

basepath = "D:/Projects/Dr. Heyl Group/data/ising/"

Temps = [1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8]
τvals = [13 , 14 , 14 , 15 , 15 , 15 , 15 , 16 , 17 , 18 , 19 ,  20]
# τvals = [3  , 4  , 8  , 11 , 20 , 185, 200, 121, 42 , 34 , 26 , 20 , 21 , 16 ]

lattice_sizes = [16, 24, 32, 40]
eqsteps = 1000
nconf = 1000
wolff=true

for N in lattice_sizes
    @time ising_getconfigdata_to_txt([N], Temps, nconf, eqsteps, τvals; store_at=basepath, wolff=wolff, ntau=10)
end