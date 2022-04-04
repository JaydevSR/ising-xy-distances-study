include("../../src/spinmc.jl")

using ScikitLearn

@sk_import neighbors: NearestNeighbors

basepath = "D:/Projects/DQCM (Dr. Heyl)/cluster_data/"
println("Using data from: $basepath")

Temps = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]

Nvals = [30, 40, 50, 60]
n_samples = 1000
xcomp=true
ycomp=true

f = Figure(resolution = (800, 600))
ax = Axis(f[1, 1], xlabel = "T", ylabel = "⟨R⋅r₁⟩ - ⟨R⟩⟨r₁⟩ ", title = "Correlation between r₁ and R with temperature")

Nmarks = Dict(10 => :circle, 20 => :rect, 30 => :diamond, 40 => :cross, 50 => :xcross, 60 => :star4)
Ncols = Dict(10 => :red, 20 => :green, 30 => :blue, 40 => :yellow, 50 => :magenta, 60 => :black)

for N in Nvals
    println("Calculating For N=$(N) ...")

    corr_arr = zeros(Float64, length(Temps))

    struc_facs_store = readdlm(basepath*"struc_facs/struc_facs_N$(N).txt", ',', Float64)

    struc_facs_store = Dict([(struc_facs_store[i, 1], struc_facs_store[i, 2:end]) for i=1:size(struc_facs_store)[1]])
    
    for stepT in 1:length(Temps)
        T = Temps[stepT]
        println("   | Temperature = $(T) ...")

        datafile = basepath*"data/xy/Size$(N)/uncorr_configs_Temp$(Temps[stepT])_N$(N).txt"
        
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
    println("   |   > Done.")
end

axislegend(ax, position = :rb)

location1 = "results/xy/r1_R_correlation_v2.pdf"

println("Saving Plots")
save(location1, f)
println("Done.")