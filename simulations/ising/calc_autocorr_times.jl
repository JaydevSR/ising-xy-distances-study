include("../../src/spinmc.jl")

store_at = "D:/Projects/Dr. Heyl Group/data/ising/"

Temps = [2.0, 2.1, 2.2, 2.24, 2.26, 2.265, 2.27, 2.275, 2.28, 2.3, 2.4, 2.5]
# Temps = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5]

lattice_sizes = [16, 24, 32]

datamat = zeros(Float64, (2, length(Temps)))
datamat[1,:] = Temps

current_loc = pwd()
cd(store_at)
for N in lattice_sizes
    println("Calculating for Size=$(N)x$(N) ...")
    location="ising_corrtimes_N$(N)_wolff.txt"
    datamat[2,:] = ising_getcorrtime(N, Temps; wolff=true, verbose=true)
    open(location, "w") do io
        writedlm(io, datamat, ',')
    end;
end
cd(current_loc)