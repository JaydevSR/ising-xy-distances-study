include("../../src/spinmc.jl")

datapath = joinpath("d:\\", "Projects", "Dr. Heyl Group", "cluster_data")
println("Using data from: $datapath")

lattice_sizes = [32, 40, 48, 56, 64]
temps = [1.800, 1.900, 2.000, 2.100, 2.200, 2.210, 2.220, 2.230, 2.240,
    2.250, 2.260, 2.270, 2.274, 2.278, 2.282, 2.286, 2.290, 2.294, 2.298,
    2.302, 2.306, 2.310, 2.314, 2.318, 2.322, 2.330, 2.340, 2.350, 2.360,
    2.370, 2.380, 2.390, 2.400, 2.500, 2.600, 2.700];

using ScikitLearn
@sk_import neighbors:NearestNeighbors
model = NearestNeighbors(n_neighbors=2, algorithm="auto", n_jobs=-1)

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
            continue
        end

        el = @elapsed begin
            configs_vecs = collect(eachcol(readdlm(datafile, ',', Int64)))
            nnbrs = fit!(model, configs_vecs)
            dists, idxs = NearestNeighbors.kneighbors(nnbrs, configs_vecs)
            fnn = dists[:, 2]

            szpath = joinpath([datapath, "Ising", "Size$L", "fnn_dists"])
            mkpath(szpath)
            open(joinpath([szpath, "ising_fnn_dists_temp$(T)_L$(L).txt"]), "w") do io
                writedlm(io, fnn, ',')
            end
        end
        println("| Process completed for (T = $(T)) in $(el) seconds.")
    end
    println("| Done.")
    println(".==================================")
end
