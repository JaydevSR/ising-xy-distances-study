include("../../src/spinmc.jl")

basepath = "D:/Projects/Dr. Heyl Group/data/ising/"
println("Using data from: $basepath")

Temps = [1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0]
Nvals = [16, 24, 32, 48]

# using NearestNeighbors

# for N in Nvals
#     @time begin
#         println("Calculating For N=$(N) ...")
#         Threads.@threads for stepT in eachindex(Temps)
#             T = Temps[stepT]
#             println("| Process strarted on thread #$(Threads.threadid()): T = $(T)")
#             datafile = basepath*"Size$(N)/ising_uncorr_configs_Temp$(Temps[stepT])_N$(N).txt"
#             if !isfile(datafile)
#                 println("| No Data Found (T=$(T)).")
#                 println("| Process ended on thread #$(Threads.threadid()): T = $(T)")
#                 continue
#             end
#             uncorrelated_spins = readdlm(datafile, ',', Float64)
#             kdtree = KDTree(uncorrelated_spins)
#             idxs, dists = knn(kdtree, uncorrelated_spins, 2, true)
#             fnn = last.(dists)
    
#             ispath(basepath*"fnn_dists/Size$N/") ? 1 : mkpath(basepath*"fnn_dists/Size$N/")
#             location=basepath*"fnn_dists/Size$N/fnn_dists_Temp$(T)_N$(N).txt"
#             open(location, "w") do io
#                 writedlm(io, fnn, ',')
#             end;
#             println("| Process completed on thread #$(Threads.threadid()): T = $(T)")
#         end        
#     end
# end


using ScikitLearn

@sk_import neighbors: NearestNeighbors

for N in Nvals
    @time begin
        println("Calculating For N=$(N) ...")
        for stepT in eachindex(Temps)
            T = Temps[stepT]
            println("| Process strarted on thread #$(Threads.threadid()): T = $(T)")
            datafile = basepath*"Size$(N)/ising_uncorr_configs_Temp$(Temps[stepT])_N$(N).txt"
            if !isfile(datafile)
                println("| No Data Found (T=$(T)).")
                println("| Process ended on thread #$(Threads.threadid()): T = $(T)")
                continue
            end
            uncorrelated_spins = reshape(readdlm(datafile, ',', Float64), N, N, :)
            model = NearestNeighbors(n_neighbors = 2, algorithm = "ball_tree")
            configs_vec = [reshape(uncorrelated_spins[:, :, i], N * N) for i = 1:size(uncorrelated_spins)[3]]
            nnbrs = fit!(model, configs_vec)
            dists, idxs = NearestNeighbors.kneighbors(nnbrs, configs_vec)
            fnn = dists[:, 2]
            ispath(basepath*"fnn_dists/Size$N/") ? 1 : mkpath(basepath*"fnn_dists/Size$N/")
            location=basepath*"fnn_dists/Size$N/fnn_dists_Temp$(T)_N$(N).txt"
            open(location, "w") do io
                writedlm(io, fnn, ',')
            end;
            println("| Process completed on thread #$(Threads.threadid()): T = $(T)")
        end
    end
end