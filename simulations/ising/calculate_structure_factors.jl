include("../../src/spinmc.jl")

basepath = "D:/Projects/Dr. Heyl Group/data/ising/"
println("Using data from: $basepath")

Temps = [1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0]


Nvals = [16, 24, 32, 48]
n_samples = 2000
for N in Nvals
    println("Calculating For N=$(N) ...")
    Threads.@threads for stepT in eachindex(Temps)
        T = Temps[stepT]
        println("| Process strarted on thread #$(Threads.threadid()): T = $(T)")
        datafile = basepath*"Size$(N)/ising_uncorr_configs_Temp$(Temps[stepT])_N$(N).txt"
        if !isfile(datafile)
            println("| No Data Found (T=$(T)).")
            println("| Process ended on thread #$(Threads.threadid()): T = $(T)")
            continue
        end
        uncorrelated_spins = reshape(readdlm(datafile, ',', Float64), N, N, :)
        struc_facs = [structure_factor(uncorrelated_spins[:, :, i], N; scaled=false) for i = 1:size(uncorrelated_spins)[3]]

        ispath(basepath*"struc_facs/Size$N/") ? 1 : mkpath(basepath*"struc_facs/Size$N/")
        location=basepath*"struc_facs/Size$N/struc_facs_Temp$(T)_N$(N).txt"
        open(location, "w") do io
            writedlm(io, struc_facs, ',')
        end;
        println("| Process completed on thread #$(Threads.threadid()): T = $(T)")
    end
end