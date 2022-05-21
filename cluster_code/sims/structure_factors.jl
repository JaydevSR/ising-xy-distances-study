include("../../src/spinmc.jl")

basepath = "D:/Projects/DQCM (Dr. Heyl)/cluster_data/data/xy/"
println("Using data from: $basepath")

Temps = [0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.975, 1.0, 1.005, 1.01, 1.015, 1.025, 1.04, 1.05, 1.075, 1.1, 1.15, 1.2, 1.25, 1.3, 1.4, 1.5]

Nvals = [60] #[10, 20, 30, 40, 50, 60]
n_samples = 1000

for N in Nvals
    println("Calculating For N=$(N) ...")
    struc_facs = zeros(Float64, (length(Temps), n_samples))
    for stepT in 1:length(Temps)
        T = Temps[stepT]
        println("   | Temperature = $(T) ...")
        
        datafile = basepath*"Size$(N)/uncorr_configs_Temp$(Temps[stepT])_N$(N).txt"

        if !isfile(datafile)
            println("   |   > No Data Found.")
            continue
        end
    
        uncorrelated_spins = reshape(readdlm(datafile, ',', Float64), N, N, :)
    
        println("   |   > Calculating Structure Factors ...")
        struc_facs[stepT, :] = [structure_factor(uncorrelated_spins[:, :, i], N; metric=xy_spindot, scaled=true) for i = 1:size(uncorrelated_spins)[3]]
        println("   |   > Done.")
    end

    ispath(basepath*"struc_facs") ? 1 : mkpath(basepath*"struc_facs")
    location=basepath*"struc_facs/struc_facs_N$(N).txt"
    open(location, "w") do io
        writedlm(io, [Temps struc_facs], ',')
    end;
end