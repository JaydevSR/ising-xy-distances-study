include("../../src/spinmc.jl")

basepath = "D:/Projects/Dr. Heyl Group/data/ising/"
println("Using data from: $basepath")

Temps = [1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0]

Nvals = [16, 24, 32, 48]
n_samples = 2000

f = Figure(resolution = (1200, 800))
ax = Axis(f[1, 1], xlabel = "T", ylabel = "C (scaled by nsites)", title = "Correlation between r₁ and Ns*R with temperature")
ax2 = Axis(f[1, 2], xlabel = "T", ylabel = "C", title = "Correlation between r₁ and R with temperature")

Nmarks = Dict(16 => :circle, 24 => :rect, 32 => :diamond, 48 => :cross, 64 => :xcross)
Ncols = Dict(16 => :red, 24 => :green, 32 => :blue, 48 => :magenta, 64 => :black)

corr_data = zeros(Float64, (length(Nvals), length(Temps)))

for stepN in 1:length(Nvals)
    N = Nvals[stepN]
    println("Calculating For N=$(N) ...")
    corr_arr = zeros(Float64, length(Temps))
    
    for stepT in eachindex(Temps)
        T = Temps[stepT]

        fnn_dists_store = basepath*"fnn_dists/Size$N/fnn_dists_Temp$(T)_N$(N).txt"
        if !isfile(fnn_dists_store)
            println("| No FNN Distance Data Found (T=$(T)).")
            println("| Process ended on thread #$(Threads.threadid()): T = $(T)")
            continue
        end
        fnn_dists = readdlm(fnn_dists_store, ',', Float64)
    
        struc_facs_store = basepath*"struc_facs/Size$N/struc_facs_Temp$(T)_N$(N).txt"
        if !isfile(struc_facs_store)
            println("| No Structure Factor Data Found (T=$(T)).")
            println("| Process ended on thread #$(Threads.threadid()): T = $(T)")
            continue
        end
        struc_facs = readdlm(struc_facs_store, ',', Float64)

        corr_arr[stepT] = mean(fnn_dists.*struc_facs) - mean(fnn_dists)*mean(struc_facs)
        println("| Process completed on thread #$(Threads.threadid()): T = $(T)")
    end

    println("|")
    println("| Adding Plot ...")
    scatter_points = Point2f.(Temps, corr_arr)
    scatterlines!(ax, scatter_points, color=Ncols[N], label="L=$(N)", marker=Nmarks[N])

    scatter_points2 = Point2f.(Temps, corr_arr/(N^2))
    scatterlines!(ax2, scatter_points2, color=Ncols[N], label="L=$(N)", marker=Nmarks[N])
    println("| Done.")
    println("|")

    corr_data[stepN, 1:length(corr_arr)] = corr_arr
end

axislegend(ax, position = :rb)
axislegend(ax2, position = :rb)

location_plot = "results/ising/r1_R_correlation_plot_v1_metro.pdf"

open(basepath*"r1_R_corr_data.txt", "w") do io
    writedlm(io, corr_data, ',')
end;

open(basepath*"r1_R_corr_Nvals.txt", "w") do io
    writedlm(io, Nvals, ',')
end;

open(basepath*"r1_R_corr_Temps.txt", "w") do io
    writedlm(io, Temps, ',')
end;

println("Saving Plots")
save(location_plot, f)
println("Done.")