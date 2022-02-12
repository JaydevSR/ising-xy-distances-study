function ising_getconfigdata!(
    file::JLD2.JLDFile,
    lattice_sizes::AbstractArray{Int64},
    Temps::AbstractArray{Float64},
    n_uncorr::Int64,
    eqsteps::Int64; 
    from_infinity::Bool=false,
    wolff::Bool=false
    )
    for N in lattice_sizes
        println("Generating configurations: Size = $(N)x$(N) ...")

        if from_infinity
            spins = rand([-1.0, 1.0], (N, N))
        else
            spins = ones(Float64, (N, N))
        end

        autocorr_times = zeros(Int64, length(Temps))
        configs_data = []
        sizehint!(configs_data, length(Temps))
        for stepT in 1:length(Temps)
            T = Temps[stepT]
            println("   | Temperature = $(T) ...")
    
            ising_equilibrate_system!(spins, T, eqsteps; wolff=wolff)

            println("   |   | Calculating correlation time ...")
            autocorr_times[stepT] = τ = ising_getcorrtime!(spins, T; wolff=wolff)
            println("   |   | Done.")

            println("   |   | Making uncorrelated measurements (τ=$(τ)) ...")
            uncorrelated_spins = ising_getuncorrconfigs!(spins, T, τ, n_uncorr; wolff=wolff)
            push!(configs_data, (T, copy(uncorrelated_spins)))
            println("   |   | Done.")
        end
    
        # Write generated data to the file
        println("\nWriting data...")
        file["$(N)x$(N)/autocorr_times"] = [(Temps[i], autocorr_times[i]) for i=1:length(Temps)]
        file["$(N)x$(N)/uncorr_configs"] = configs_data
    
        println("Done.")
        println("------------------------------------------------\n")
    end
end

function ising_getcorrtime!(spins::Matrix, T, msteps=6000; wolff=false)
    N = size(spins)[1]
    mags = zeros(Float64, msteps)
    E0 = ising_total_energy(spins)
    mags[1] = ising_total_magnetization(spins)
    if wolff
        P_add = isingwolff_Padd(T)
        for i in 1:msteps-1
            isingwolff_step!(spins, P_add)
            mags[i+1] = ising_total_magnetization(spins)
        end
    else
        for i in 1:msteps-1
            E0, mags[i+1] = isingmetro_step!(spins, T, E0, mags[i])
        end
    end
    
    @. mags /= N^2
    corrfn = autocorrelation_fn(mags)

    # Measure the integrated correlation time in a small window
    τ = corrfn[1:200] |> sum |> ceil
    return convert(Int64, τ)
end

function ising_getuncorrconfigs!(spins::Matrix, T, τ, n_uncorr; wolff=false)
    N = size(spins)[1]
    twice_τ = 2*τ
    nsteps = twice_τ*n_uncorr
    uncorrelated_spins = zeros(Float64, (N, N, n_uncorr))
    # uncorrelated measurements
    E0, M0 = ising_total_energy(spins), ising_total_magnetization(spins)

    if wolff
        P_add = isingwolff_Padd(T)
        for j=1:nsteps
            isingwolff_step!(spins, P_add)
            if j%twice_τ == 0
                uncorrelated_spins[:, :, j÷twice_τ] = spins
            end
        end
    else
        for j in 1:nsteps
            E0, M0 = isingmetro_step!(spins, T, E0, M0)
            if j%twice_τ == 0
                uncorrelated_spins[:, :, j÷twice_τ] = spins
            end
        end
    end
    return uncorrelated_spins
end

function ising_equilibrate_system!(spins::Matrix, T, eqsteps; wolff=false)
    E0, M0 = ising_total_energy(spins), ising_total_magnetization(spins)
    if wolff
        for i in 1:eqsteps
            P_add = isingwolff_Padd(T)
            isingwolff_step!(spins, P_add)
        end
    else
        for i in 1:eqsteps
            E0, M0 = isingmetro_step!(spins, T, E0, M0)
        end
    end
end