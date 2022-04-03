include("../../src/spinmc.jl")

using ScikitLearn

@sk_import neighbors: NearestNeighbors

datapath = "results/xy/xy_config_data.jld2"
println("Using data from: $datapath")
file = jldopen(datapath)

Temps = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]

Nvals = [10, 20, 30, 40]
n_samples = 500
xcomp=true
ycomp=true

f = Figure(resolution = (800, 600))
f2 = Figure(resolution = (800, 600))
ax = Axis(f[1, 1], xlabel = "⟨r₁⟩", ylabel = "Γ", title = "Spin-stiffness vs mean value of r₁ ")
ax2 = Axis(f2[1, 1], xlabel = "T", ylabel = "Γ", title = "Spin-stiffness vs Temperature")

Nmarks = Dict(10 => :circle, 20 => :rect, 30 => :diamond, 40 => :cross)

for N in Nvals
    println("Calculating For N=$(N) ...")
    data = Dict(file["$(N)x$(N)/uncorr_configs"])

    mean_r1 = zeros(Float64, length(Temps))
    spin_stiffness = zeros(Float64, length(Temps))
    
    for stepT in 1:length(Temps)
        T = Temps[stepT]
        println("   | Temperature = $(T) ...")
        
        if haskey(data, T)
            uncorrelated_spins = data[T]
        else
            println("   |   > No Data Found.")
            continue
        end
    
        uncorrelated_spins = data[T]
    
        println("   |   > Calculating Distances ...")
        # fit the kNN algorithm to uncorrelated spin configs
        model = NearestNeighbors(n_neighbors = 2, algorithm = "ball_tree")
        configs_vec = xy_prepare_vector(uncorrelated_spins; xcomp=xcomp, ycomp=ycomp)
        nnbrs = fit!(model, configs_vec)  # Ignore warning 
    
        # calculate distances
        dists, idxs = NearestNeighbors.kneighbors(nnbrs, configs_vec)
        println("   |   > Done.")
        mean_r1[stepT] = mean(dists[:, 2])
    
        println("   |   > Calculating Spin-stiffness ...")
        configs = [uncorrelated_spins[:, :, i] for i in 1:n_samples]
        spin_stiffness[stepT] = xy_spin_stiffness_value(configs, T, N)
        println("   |   > Done.")
    end

    println("   |   > Adding Plot ...")
    scatter_points = Point2f.(mean_r1, spin_stiffness)
    scatterlines!(ax, scatter_points, color=1:length(Temps), label="L=$(N)", marker=Nmarks[N])
    scatterlines!(ax2, Temps, spin_stiffness, label="L=$(N)", color=cgrad(:default, 100, categorical = true)[rand(1:100)], marker=Nmarks[N])
    println("   |   > Done.")
end

Colorbar(
    f[1,2],
    # colormap = cgrad(:matter, rev=true),
    limits = (Temps[1], Temps[end]),
    ticks = Temps[1:end], tickalign=1,
    label="Temperature"
)
axislegend(ax)
axislegend(ax2)

location1 = "results/xy/mean_r1_spin_stiffness.pdf"
location2 = "results/xy/spin_stiffness_with_temperature.pdf"

println("Saving Plots")
save(location1, f)
save(location2, f2)
println("Done.")

close(file)