include("../../src/spinmc.jl")

datapath = joinpath("d:\\", "Projects", "Dr. Heyl Group", "cluster_data")
println("Using data from: $datapath")

lattice_sizes = [16, 32, 48, 64];
temps = [1.80, 1.90, 2.00, 2.10, 2.20, 2.21, 2.22, 2.23, 2.24, 2.25, 2.26, 2.27, 2.28, 2.29, 2.30, 2.31, 2.32, 2.33, 2.34, 2.35, 2.36, 2.37, 2.38, 2.39, 2.40, 2.50, 2.60, 2.70];

using ScikitLearn
@sk_import neighbors: NearestNeighbors
model = NearestNeighbors(n_neighbors = 2, algorithm = "ball_tree")

for L in lattice_sizes
    println(".==================================")
    println("| Lattice Size: $(L) x $(L)        ")
    println(".==================================")
    println("|  ")
    for stepT in eachindex(temps)
        T = temps[stepT]
        println("| Process strarted for (T = $(T)).")
        datafile = joinpath([datapath, "Ising", "Size$L", "uncorr_configs", "ising_uncorr_configs_temp$(T)_L$(L).txt"])
        if !isfile(datafile)
            println("| No Data Found (T=$(T)).")
            println("| Process ended on thread #$(Threads.threadid()): T = $(T)")
            continue
        end

        el = @elapsed begin
            configs_vecs = collect(eachcol(readdlm(datafile, ',', Float64)))
            nnbrs = fit!(model, configs_vecs)
            dists, idxs = NearestNeighbors.kneighbors(nnbrs, configs_vecs)
            fnn = dists[:, 2]
    
            szpath = joinpath([datapath, "Ising", "Size$L", "fnn_dists"])
            mkpath(szpath)
            open(joinpath([szpath, "ising_fnn_dists_temp$(T)_L$(L).txt"]), "w") do io
                writedlm(io, fnn, ',')
            end; 
        end
        println("| Process completed for (T = $(T)) in $(el) seconds.")
    end
    println("| Done.")
    println(".==================================")
end
