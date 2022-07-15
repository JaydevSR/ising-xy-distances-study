using DelimitedFiles
using CairoMakie
using GLM
using DataFrames

data_path = "D:/Projects/Dr. Heyl Group/data/ising/"

Temps = readdlm(data_path*"temperature_values.txt", ',', Float64)[1:end,]
Nvals = readdlm(data_path*"lattice_sizes.txt", ',', Float64)[1:end,]
corr_data = readdlm(data_path*"r1_R_correlation_data.txt", ',', Float64)

T_c = 2 / log1p(sqrt(2))

ln_arr = 1 ./ Nvals
T_star_arr = zeros(Float64, length(Nvals))

for i=1:length(Nvals)
    T_star_arr[i] = Temps[argmax(abs.(corr_data[i, :]))]
end

T_star_arr = T_star_arr #.- T_c

data_all = DataFrame(X=ln_arr[end-2:2:end], Y=T_star_arr[end-2:2:end])

ols1 = lm(@formula(Y~X), data_all)

y_reg = coef(ols1)[2].*append!([0.0], ln_arr) .+ coef(ols1)[1]

location_plot = "results/ising/"
f = Figure(resolution = (800, 600))
ax = Axis(f[1,1], xlabel = "1/ln(L)^2", ylabel = "T*", title="Scaling of T* (T_c = $T_c). Fit: y = $(coef(ols1)[2])x + $(coef(ols1)[1])")

scatter_points = Point2f.(ln_arr, T_star_arr)
scatter!(ax, scatter_points, marker=:diamond, color=:red, markersize=14)
CairoMakie.xlims!(ax, (0, 0.07))
CairoMakie.ylims!(ax, (2.2, 2.35))



lines!(ax, Point2f.(append!([0.0], ln_arr), y_reg), label="liner fit", color=:black)

# save(location_plot*"T_star_scaling_plot_regression_v1.pdf", f)
display(f)
