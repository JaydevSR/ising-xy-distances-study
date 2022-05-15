using DelimitedFiles
# using Plots
using CairoMakie

data_path = "D:/Projects/DQCM (Dr. Heyl)/montecarlo-spin-configurations-and-distances-study/results/xy/"

Temps = reshape(readdlm(data_path*"corr_Temps.txt", ',', Float64), (:,))
Nvals = reshape(readdlm(data_path*"corr_Nvals.txt", ',', Float64), (:,))
corr_data = readdlm(data_path*"corr_data.txt", ',', Float64)

# Plot T* vs 1/(ln(L))^2

ln_arr = 1 ./ ((log.(Nvals)).^2)
T_star_arr = zeros(Float64, length(Nvals))

for i=1:length(Nvals)
    T_star_arr[i] = Temps[argmax(abs.(corr_data[i, :]))]
end

location_plot = "results/xy/"
f = Figure(resolution = (800, 600))
ax = Axis(f[1,1], xlabel = "1/ln(L)^2", ylabel = "T*")

println(ln_arr)
println(T_star_arr)

scatter_points = Point2f.(ln_arr, T_star_arr)
scatter!(ax, scatter_points, marker=:diamond, color=:red, markersize=14)
xlims!(ax, (0, 0.125))
ylims!(ax, (0.85, 1.15))
save(location_plot*"T_star_scaling_plot.pdf", f)