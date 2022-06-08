using DelimitedFiles
using CairoMakie
using LaTeXStrings

data_path = "D:/Projects/DQCM (Dr. Heyl)/montecarlo-spin-configurations-and-distances-study/results/xy/"

Temps = reshape(readdlm(data_path*"corr_Temps.txt", ',', Float64), (:,))
Lvals = reshape(readdlm(data_path*"corr_Nvals.txt", ',', Float64), (:,))
corr_data = readdlm(data_path*"corr_data.txt", ',', Float64)

# BKT universality class

T_KT = 0.8933
b = 1.7

# ξ ~ exp(b*t^-0.5)
# t = T - T_KT

p = 2.26

start_idx = 1
# T > T_KT
for i in 1:length(Temps)
    global start_idx
    if Temps[i] > T_KT
        start_idx = i
        break
    end
end

t_vals = Temps[start_idx:end] .- T_KT

location_plot = "results/xy/"
f = Figure(resolution = (800, 600))
ax = Axis(f[1,1], xlabel = L"L\exp(-bt^{-0.5})", ylabel = L"L^{-p}C_{sft, r_1}", title=L"data collpse of correlation data: $T_{KT} = %$T_KT$, $b = %$b$, $p = %$p$")

Lmarks = Dict(10 => :circle, 20 => :rect, 30 => :diamond, 40 => :cross, 50 => :xcross, 60 => :star4)
Lcols = Dict(10 => :red, 20 => :green, 30 => :blue, 40 => :yellow, 50 => :magenta, 60 => :black)

for i in 1:length(Lvals)
    L = Lvals[i]
    scaling_variable = L .* exp.(-b ./ (abs.(t_vals).^0.5))
    scaling_fn = corr_data[i, start_idx:end] ./ L^p
    scatter!(ax, Point2f.(scaling_variable, scaling_fn), marker=Lmarks[L], color=Lcols[L], label="L=$(L)")
end

axislegend(ax, position=:rb)

save(location_plot*"corr_data_collapse_plot_v1.pdf", f)