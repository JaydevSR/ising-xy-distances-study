include("../../src/spinmc.jl")

basepath = "D:/Projects/Dr. Heyl Group/data/ising/"
println("Using data from: $basepath")

Temps = [1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0]

Nvals = [16, 24, 32, 48]
nconf = 1000

f=Figure()
ax=Axis(f[1,1], xlabel="temperature", ylabel="mean magnetisation")
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
        uncorrelated_spins = reshape(readdlm(datafile, ',', Float64), N, N, :)[:, :, 1:nconf]
        
        for i=1:nconf
            mag[stepT] += abs(ising_total_magnetization(uncorrelated_spins[:, :, i]))
        end
        mag[stepT] /= N*N*nconf
        println("| Process completed on thread #$(Threads.threadid()): T = $(T)")
    end
    scatterlines!(Temps, mag, label="L=$(N)")
end

display(f)