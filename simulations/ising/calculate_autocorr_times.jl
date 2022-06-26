include("../../src/spinmc.jl")

wolff = true
msteps = 50000
wsteps = 1000

basepath = "D:/Projects/Dr. Heyl Group/data/"

Temps = readdlm(joinpath(basepath, "ising",  "temperature_values.txt"), ',', Float64)[:, 1]

lattice_sizes = readdlm(joinpath(basepath, "ising", "lattice_sizes.txt"), ',', Int64)[:, 1]

for N in lattice_sizes
    println("Calculating for Size=$(N)x$(N) ...")
    location=joinpath(basepath, "ising", "autocorr_times", "ising_autocorr_times_size$(N)_$(wolff ? "wolff" : "metro").txt")
    autocorr_times = ising_getcorrtime(N, Temps, msteps, wsteps; wolff=wolff, verbose=true)
    open(location, "w") do io
        writedlm(io, autocorr_times, ',')
    end;
end