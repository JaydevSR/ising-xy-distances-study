include("../../src/spinmc.jl")

store_at = "D:/Projects/Dr. Heyl Group/data/ising/"

Temps = [2.0, 2.1, 2.2, 2.24, 2.26, 2.265, 2.27, 2.275, 2.28, 2.3, 2.4, 2.5]

lattice_sizes = [16, 24, 32]
eqsteps = 10000
n_uncorr = 10000

current_loc = pwd()
cd(store_at)
for N in lattice_sizes
    println("Calculating for Size=$(N)x$(N) ...")
    location="ising_corrtimes_N$(N)_wolff.txt"
    corr_times = ising_getcorrtime(N, Temps; wolff=true, verbose=true)
    open(location, "w") do io
        writedlm(io, corr_times, ',')
    end;
end
cd(current_loc)