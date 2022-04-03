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

function xy_getconfigdata!(
    file::JLD2.JLDFile,
    lattice_sizes::AbstractArray{Int64},
    Temps::AbstractArray{Float64},
    n_uncorr::Int64,
    eqsteps::Int64; 
    calculate_autocorr_times::Bool=false,
    autocorr_times::AbstractArray{Int64}=ones(Int64, length(Temps))
    )
    for N in lattice_sizes
        println("Generating configurations: Size = $(N)x$(N) ...")

        spins = rand(Float64, (N, N))

        configs_data = []
        sizehint!(configs_data, length(Temps))
        for stepT in 1:length(Temps)
            T = Temps[stepT]
            println("   | Temperature = $(T) ...")
    
            xy_equilibrate_system!(spins, T, eqsteps)

            if calculate_autocorr_times
                println("   |   | Calculating correlation time ...")
                autocorr_times[stepT] = xy_getcorrtime!(spins, T)
                println("   |   | Done.")
            end

            τ = autocorr_times[stepT]
            println("   |   | Making uncorrelated measurements (τ=$(τ)) ...")
            uncorrelated_spins = xy_getuncorrconfigs!(spins, T, τ, n_uncorr)
            push!(configs_data, (T, copy(uncorrelated_spins)))
            println("   |   | Done.")
        end
    
        # Write generated data to the file
        println("\nWriting data...")
        # file["$(N)x$(N)/autocorr_times"] = [(Temps[i], autocorr_times[i]) for i=1:length(Temps)]
        file["$(N)x$(N)/uncorr_configs"] = configs_data
    
        println("Done.")
        println("------------------------------------------------\n")
    end
end

function xy_getconfigdata_to_txt(
    lattice_sizes::AbstractArray{Int64},
    Temps::AbstractArray{Float64},
    n_uncorr::Int64,
    eqsteps::Int64; 
    calculate_autocorr_times::Bool=false,
    autocorr_times::AbstractArray{Int64}=ones(Int64, length(Temps)),
    store_at::AbstractString=""
    )
    current_loc = pwd()
    cd(store_at)
    for N in lattice_sizes
        println("Generating configurations: Size = $(N)x$(N) ...")
        spins = rand(Float64, (N, N))
        for stepT in 1:length(Temps)
            isdir("Size$N") ? 1 : mkdir("Size$N")
            cd("Size$N")
            location="uncorr_configs_Temp$(Temps[stepT])_N$(N).txt"

            file = open(location, "w")
            T = Temps[stepT]
            println("   | Temperature = $(T) ...")
    
            xy_equilibrate_system!(spins, T, eqsteps)

            if calculate_autocorr_times
                println("   |   | Calculating correlation time ...")
                autocorr_times[stepT] = xy_getcorrtime!(spins, T)
                println("   |   | Done.")
            end

            τ = autocorr_times[stepT]
            println("   |   | Making uncorrelated measurements (τ=$(τ)) ...")
            uncorrelated_spins = xy_getuncorrconfigs!(spins, T, τ, n_uncorr)
            
            open(location, "w") do io
                writedlm(io, reshape(uncorrelated_spins, (N*N, n_uncorr)), ',')
            end;
            cd("..")
            println("   |   | Done.")
        end
        cd(current_loc)
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

function xy_getcorrtime!(spins::Matrix, T, msteps=5000)
    N = size(spins)[1]
    ergs = zeros(Float64, msteps)
    for i in 1:msteps-1
        xywolff_step!(spins, N, T)
        ergs[i+1] = xy_total_energy(spins, N)
    end
    
    @. ergs /= N^2
    corrfn = autocorrelation_fn(ergs)

    # Measure the integrated correlation time in a small window
    τ = corrfn[1:200] |> sum |> ceil
    return convert(Int64, τ)
end

function xy_getuncorrconfigs!(spins::Matrix, T, τ, n_uncorr)
    N = size(spins)[1]
    twice_τ = 2*τ
    nsteps = twice_τ*n_uncorr
    uncorrelated_spins = zeros(Float64, (N, N, n_uncorr))
    for j=1:nsteps
        xywolff_step!(spins, N, T)
        if j%twice_τ == 0
            uncorrelated_spins[:, :, j÷twice_τ] = spins
        end
    end
    return uncorrelated_spins
end

function xy_equilibrate_system!(spins::Matrix, T, eqsteps)
    N = size(spins)[1]
    for i in 1:eqsteps
        xywolff_step!(spins, N, T)
    end
end

function xy_prepare_vector(config_array::AbstractArray; xcomp=true, ycomp=true)
    nconf = size(config_array)[3]
    N = size(config_array)[1]
    configs_vec = [reshape(config_array[:, :, i], N * N) for i = 1:nconf]
    emptyness = []
    if xcomp
        x_vec = map(x -> cos2pi.(x), configs_vec)
    else
        x_vec = fill(emptyness, nconf)
    end

    if ycomp
        y_vec = map(y -> sin2pi.(y), configs_vec)
    else
        y_vec = fill(emptyness, nconf)
    end

    spins_vector = []
    sizehint!(spins_vector, nconf)
    for i=1:nconf
        push!(spins_vector, [x_vec[i]..., y_vec[i]...])
    end
    return spins_vector
end

function get_array_from_file(file, N)
    arr = readdlm("file", ",", Float64)
    return reshape(arr, N, N, :)
end