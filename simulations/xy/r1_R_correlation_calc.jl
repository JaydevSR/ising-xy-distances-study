include("../../src/spinmc.jl")

using ScikitLearn

@sk_import neighbors: NearestNeighbors

basepath = "D:/Projects/DQCM (Dr. Heyl)/cluster_data/data/xy/"
println("Using data from: $basepath")

Temps = [0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.975, 1.0, 1.01, 1.025, 1.04, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5]

Nvals = [10, 20, 30, 40, 50, 60]
n_samples = 1000
xcomp=true
ycomp=true

f = Figure(resolution = (1200, 800))
ax = Axis(f[1, 1], xlabel = "T", ylabel = "C (scaled by nsites)", title = "Correlation between r₁ and Ns*R with temperature")
ax2 = Axis(f[1, 2], xlabel = "T", ylabel = "C", title = "Correlation between r₁ and R with temperature")

Nmarks = Dict(10 => :circle, 20 => :rect, 30 => :diamond, 40 => :cross, 50 => :xcross, 60 => :star4)
Ncols = Dict(10 => :red, 20 => :green, 30 => :blue, 40 => :yellow, 50 => :magenta, 60 => :black)

corr_data = zeros(Float64, (length(Nvals), length(Temps)))

for stepN in 1:length(Nvals)
    N = Nvals[stepN]
    println("Calculating For N=$(N) ...")

    corr_arr = zeros(Float64, length(Temps))

    struc_facs_store = readdlm(basepath*"struc_facs/struc_facs_N$(N).txt", ',', Float64)

    struc_facs_store = Dict([(struc_facs_store[i, 1], struc_facs_store[i, 2:end]) for i=1:size(struc_facs_store)[1]])
    
    for stepT in 1:length(Temps)
        T = Temps[stepT]
        println("   | Temperature = $(T) ...")

        datafile = basepath*"Size$(N)/uncorr_configs_Temp$(Temps[stepT])_N$(N).txt"
        
        if !isfile(datafile)
            println("   |   > No Data Found.")
            continue
        end
    
        uncorrelated_spins = reshape(readdlm(datafile, ',', Float64), N, N, :)
    
        println("   |   > Calculating Distances ...")
        # fit the kNN algorithm to uncorrelated spin configs
        model = NearestNeighbors(n_neighbors = 2, algorithm = "ball_tree")
        configs_vec = xy_prepare_vector(uncorrelated_spins; xcomp=xcomp, ycomp=ycomp)
        nnbrs = fit!(model, configs_vec)  # Ignore warning 
    
        # calculate distances
        dists, idxs = NearestNeighbors.kneighbors(nnbrs, configs_vec)
        println("   |   > Done.")
    
        println("   |   > Fetching Structure Factors ...")
        if !haskey(struc_facs_store, Temps[stepT])
            println("   |   |   > No Temperature Data Found.")
            continue
        end

        struc_facs = struc_facs_store[Temps[stepT]]
        println("   |   > Done.")

        corr_arr[stepT] = mean(dists.*struc_facs) - mean(dists)*mean(struc_facs)
    end

    println("   |   > Adding Plot ...")
    scatter_points = Point2f.(Temps, corr_arr)
    scatterlines!(ax, scatter_points, color=Ncols[N], label="L=$(N)", marker=Nmarks[N])

    scatter_points2 = Point2f.(Temps, corr_arr/(N^2))
    scatterlines!(ax2, scatter_points2, color=Ncols[N], label="L=$(N)", marker=Nmarks[N])
    println("   |   > Done.")

    corr_data[stepN, 1:length(corr_arr)] = corr_arr
end

axislegend(ax, position = :rb)
axislegend(ax2, position = :rb)

location_plot = "results/xy/r1_R_correlation_plot.pdf"
location_data = "results/xy/"

open(location_data*"corr_data.txt", "w") do io
    writedlm(io, corr_data, ',')
end;

open(location_data*"corr_Nvals.txt", "w") do io
    writedlm(io, Nvals, ',')
end;

open(location_data*"corr_Temps.txt", "w") do io
    writedlm(io, Temps, ',')
end;

println("Saving Plots")
save(location_plot, f)
println("Done.")