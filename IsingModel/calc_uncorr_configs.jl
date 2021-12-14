using JLD2

include("isingmetro.jl")

isingdata = jldopen("data/isingdata.jld2", "a+")

# Simulation Parameters
lattice_sizes = [10, 20, 30, 40, 50]

eqsteps = 2000  # Steps for equilibration
N_uncorr = 50  # Number of uncorrelated measurements

for N in lattice_sizes
    println("Determining uncorrelated configurations for $(N)x$(N) lattice ...")
    spins = spins = rand([-1.0, 1.0], (N, N))
    uncorrelated_spins = zeros(N, N, N_uncorr);

    corr_time_data = isingdata["$(N)x$(N)/corr_times"]
    temps, corr_time = first.(corr_time_data), last.(corr_time_data)

    measure_times = 2 .* convert.(Int64, ceil.(corr_time))  # TODO: add this part to other file

    temp_data = []

    for stepT in 1:length(temps)
        T = temps[stepT]
        println("   | For Temperature = $(T) ...")

        # equilibration
        E0, M0 = total_energy(spins), total_magnetization(spins)
        for i in 1:eqsteps
            E0, M0 = ising_metropolis_sweep!(spins, T, E0, M0)
        end
    
        twice_τ = measure_times[stepT]
        nsteps = twice_τ * N_uncorr
        # uncorrelated measurements
        for j in 1:nsteps
            E0, M0 = ising_metropolis_sweep!(spins, T, E0, M0)
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