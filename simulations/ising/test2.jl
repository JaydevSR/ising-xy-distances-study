include("../../src/spinmc.jl")

basepath = "D:/Projects/Dr. Heyl Group/data/ising__wolff/"
println("Using data from: $basepath")

# Temps = [2.0, 2.1, 2.2, 2.24, 2.26, 2.265, 2.27, 2.275, 2.28, 2.3, 2.4]
Temps = [1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0]


Nvals = [16]
n_samples = 2000

f=Figure()
ax=Axis(f[1,1])
for N in Nvals
    println("Calculating For N=$(N) ...")
    mag = zeros(Float64, length(Temps))
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
        mag[stepT] = mean([ising_total_magnetization(uncorrelated_spins[:, :, i]) for i=1:size(uncorrelated_spins)[3]])
        println(mag)
        println("| Process completed on thread #$(Threads.threadid()): T = $(T)")
    end
    scatterlines!(ax, Point2f.(Temps, mag))
end

display(f)