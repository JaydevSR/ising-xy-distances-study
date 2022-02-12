using CairoMakie
using JLD2

data = jldopen("results/ising/ising_configdata.jld2", "r")
Temps = first.(data["10x10/mean_r1"])
mean_r1_10 = last.(data["10x10/mean_r1"])
mean_r1_20 = last.(data["20x20/mean_r1"])
f = Figure(resolution = (800, 400), font_size = 12)

ax1 = Axis(f[1,1], xlabel="Temperature", ylabel="⟨r₁⟩", title="Mean value of first NN distance with Temperature")

scatterlines!(ax1, Temps, mean_r1_10, label="Size=10x10", markercolor=:blue)
scatterlines!(ax1, Temps, mean_r1_20, label="Size=20x20", markercolor=:purple)
axislegend()

# display(f)
save("results/ising/mean_r1_with_T.png", f)

close(data)