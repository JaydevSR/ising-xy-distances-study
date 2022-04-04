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
        println("Done.")
        println("------------------------------------------------\n")
    end
    cd(current_loc)
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

function get_configs_from_file(file, N)
    arr = readdlm("file", ",", Float64)
    return reshape(arr, N, N, :)
end