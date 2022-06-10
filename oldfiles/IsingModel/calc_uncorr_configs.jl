using JLD2

# Simulation Code
include("../../src/spinmc.jl")

# Open file in read + append
isingdata = jldopen("oldfiles/IsingModel/data/isingdata.jld2", "w")

# Simulation Parameters
lattice_sizes = [10, 20]
eqsteps = 2000  # Steps for equilibration
N_uncorr = 100  # Number of uncorrelated measurements

for N in lattice_sizes
    println("Determining uncorrelated configurations for $(N)x$(N) lattice ...")

    # # Initial spins at T = ∞
    spins = spins = rand([-1.0, 1.0], (N, N))
    uncorrelated_spins = zeros(N, N, N_uncorr)

    # corr_time_data = isingdata["$(N)x$(N)/corr_times"]
    temps, corr_times = collect(1.4:0.2:3.6), fill(4, length(collect(1.4:0.2:3.6)))

    # measure configs at 2τ
    measure_times = 2 .* corr_times

    temp_data = []
    for stepT in 1:length(temps)
        T = temps[stepT]
        println("   | For Temperature = $(T) ...")

        # equilibration
        E0, M0 = ising_total_energy(spins), ising_total_magnetization(spins)
        for i in 1:eqsteps
            #E0, M0 = isingmetro_step!(spins, T, E0, M0)
            isingwolff_step!(spins, isingwolff_Padd(T))
        end
    
        twice_τ = measure_times[stepT]
        nsteps = twice_τ * N_uncorr
        
        # uncorrelated measurements
        for j in 1:nsteps
            #E0, M0 = isingmetro_step!(spins, T, E0, M0)
            isingwolff_step!(spins, isingwolff_Padd(T))
            if j%twice_τ == 0
                uncorrelated_spins[:, :, j÷twice_τ] = spins
            end
        end
        push!(temp_data, (T, copy(uncorrelated_spins)))
        println("   | Done.")
    end

    # Write generated data to the file
    println("Writing data...")
    isingdata["$(N)x$(N)/uncorr_measurements"] = temp_data

    println("Done.")
    println("----------------------------------------------------------------\n")
end

close(isingdata)