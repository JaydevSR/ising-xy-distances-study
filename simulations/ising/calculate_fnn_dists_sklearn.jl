include("../../src/spinmc.jl")

basepath = "D:/Projects/Dr. Heyl Group/data/"
println("Using data from: $basepath")

Temps = readdlm(joinpath(basepath, "temperature_values.txt"), ',', Float64)[:, 1]

lattice_sizes = readdlm(joinpath(basepath, "lattice_sizes.txt"), ',', Int64)[:, 1]

using ScikitLearn
@sk_import neighbors: NearestNeighbors
model = NearestNeighbors(n_neighbors = 2, algorithm = "ball_tree")

for N in lattice_sizes
    println(".==================================")
    println("| Lattice Size: $(N) x $(N)        ")
    println(".==================================")
    println("|  ")
    for stepT in eachindex(Temps)
        T = Temps[stepT]
        println("| Process strarted on thread #$(Threads.threadid()) (T = $(T)).")
        datafile = joinpath([basepath, "ising", "uncorr_configs", "Size$N", "ising_uncorr_configs_temp$(T)_size$(N).txt"])
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
    
            szpath = joinpath([basepath, "ising", "fnn_dists", "Size$N"])
            ispath(szpath) ? 1 : mkpath(szpath)
            open(joinpath([szpath, "fnn_dists_temp$(T)_size$(N).txt"]), "w") do io
                writedlm(io, fnn, ',')
            end; 
        end
        println("| Process completed on thread #$(Threads.threadid()) (T = $(T)) in $(el) seconds.")
    end
    println("| Done.")
    println(".==================================")
end