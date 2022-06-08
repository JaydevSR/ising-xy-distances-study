using DelimitedFiles
using CairoMakie
using GLM
using DataFrames

data_path = "D:/Projects/Dr. Heyl Group/data/ising/"

Temps = reshape(readdlm(data_path*"corr_Temps.txt", ',', Float64), (:,))
Nvals = reshape(readdlm(data_path*"corr_Nvals.txt", ',', Float64), (:,))
corr_data = readdlm(data_path*"corr_data.txt", ',', Float64)

T_c = 2 / log1p(sqrt(2))
# Plot T* vs 1/(ln(L))^2

ln_arr = 1 ./ ((log.(Nvals)).^2)
T_star_arr = zeros(Float64, length(Nvals))

for i=1:length(Nvals)
    T_star_arr[i] = Temps[argmax(abs.(corr_data[i, :]))]
end

T_star_arr = T_star_arr #.- T_c

data_all = DataFrame(X=ln_arr, Y=T_star_arr)

ols1 = lm(@formula(Y~X), data_all)

y_reg = coef(ols1)[2].*append!([0.0], ln_arr, [0.2]) .+ coef(ols1)[1]

location_plot = "results/ising/"
f = Figure(resolution = (800, 600))
ax = Axis(f[1,1], xlabel = "1/ln(L)^2", ylabel = "T*", title="Scaling of T* (T_c = $T_c). Fit: y = $(coef(ols1)[2])x + $(coef(ols1)[1])")

scatter_points = Point2f.(ln_arr, T_star_arr)
scatter!(ax, scatter_points, marker=:diamond, color=:red, markersize=14)
xlims!(ax, (0, 0.2))
ylims!(ax, (0.0, 0.4))



lines!(ax, Point2f.(append!([0.0], ln_arr, [0.2]), y_reg_1), label="liner fit", color=:black)
lines!(ax, Point2f.(append!([0.0], ln_arr, [0.2]), y_reg_2), label="liner fit (removing outlier)", color=:blue)
axislegend(ax, position=:lt)

save(location_plot*"T_star_scaling_plot_regression_v1.pdf", f)
# display(f)
