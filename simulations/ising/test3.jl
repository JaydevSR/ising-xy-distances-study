include("../../src/spinmc.jl")

N=8
T=2.2
msteps=100000
wsteps=1000
eqsteps=2000
wolff=false
spins = fill(1.0, (N, N))
mags = zeros(Float64, msteps)

if wolff
    P_add = isingwolff_Padd(T)
    for i in 1:eqsteps
        isingwolff_step!(spins, P_add)
    end
    mags[1] = ising_total_magnetization(spins)
    for i in 1:msteps-1
        isingwolff_step!(spins, P_add)
        mags[i+1] = ising_total_magnetization(spins)
    end
else
    local E0 = ising_total_energy(spins)
    mags[1] = ising_total_magnetization(spins)
    for i in 1:eqsteps
        E0, mags[1] = isingmetro_step!(spins, T, E0, mags[i])
    end
    for i in 1:msteps-1
        E0, mags[i+1] = isingmetro_step!(spins, T, E0, mags[i])
    end
end

@. mags /= N^2
corrfn = autocorrelation_fn(abs.(mags))

f = Figure()
ax = Axis(f[1,1])
s = scatter!(ax, corrfn[1:1000])
save("corrfn.png", f)