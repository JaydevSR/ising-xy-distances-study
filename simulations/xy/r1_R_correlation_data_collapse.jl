using DelimitedFiles
using CairoMakie

data_path = "D:/Projects/DQCM (Dr. Heyl)/montecarlo-spin-configurations-and-distances-study/results/xy/"

Temps = reshape(readdlm(data_path*"corr_Temps.txt", ',', Float64), (:,))
Nvals = reshape(readdlm(data_path*"corr_Nvals.txt", ',', Float64), (:,))
corr_data = readdlm(data_path*"corr_data.txt", ',', Float64)

# BKT universality class

T_c = 0.8933