using JLD2

include("../../src/spinmc.jl")

# open/overwrite file
file = jldopen("oldfiles/IsingModel/data/isingdata.jld2", "w")

# Simulation Parameters
lattice_sizes = [10, 20]
Temps = collect(1.4:0.2:3.6)
corr_times = zeros(Int64, length(Temps))

eqsteps = 2000  # Steps for equilibration
msteps = 6000  # Steps for measurements

for N in lattice_sizes
    println("Calculating correlation times for $(N)x$(N) ...")

    # Initial spins at T = ∞
    spins = rand([-1.0, 1.0], (N, N))

    for stepT in 1:length(Temps)
        T = Temps[stepT]
        println("   | Calculating for Temperature = $(T) ...")

        # Equilibration
        E0, M0 = ising_total_energy(spins), ising_total_magnetization(spins)
        for i in 1:eqsteps
            E0, M0 = isingmetro_step!(spins, T, E0, M0)
        end
        
        mags = zeros(Float64, msteps)
        mags[1] = M0
        for i in 1:msteps-1
            E0, mags[i+1] = isingmetro_step!(spins, T, E0, mags[i])
        end
        
        @. mags /= N^2
        corrfn = autocorrelation_fn(mags)

        # Measure the integrated correlation time in a small window
        corr_times[stepT] = ceil(Int64, sum(corrfn[1:200]))

        println("   | Done (τ=$(corr_times[stepT])).")
    end

    # Write generated data to the file
    println("Writing data...")
    # convert correlation times to integer values
    file["$(N)x$(N)/corr_times"] = [
        (Temps[i], convert(Int64, ceil(corr_times[i]))) 
        for i=1:length(Temps)
    ]

    println("Done.")
    println("------------------------------------------------\n")
end

close(file)