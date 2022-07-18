function ising_getconfigdata_to_txt(
    lattice_sizes::AbstractArray{Int64},
    Temps::AbstractArray{Float64},
    n_uncorr::Int64,
    eqsteps::Int64,
    autocorr_times::AbstractArray{Int64}=ones(Int64, length(Temps));
    store_at::AbstractString="",
    from_infinity::Bool=false,
    verbose=true,
    ntau=2,
    mode="w"
    )
    for N in lattice_sizes
        verbose && println(".==================================")
        verbose && println("| Lattice Size: $(N) x $(N)        ")
        verbose && println(".==================================")
        verbose && println("|  ")
        szpath = joinpath([store_at, "ising", "uncorr_configs", "Size$N"])
        ispath(szpath) ? 1 : mkpath(szpath)
        @sync for stepT in 1:length(Temps)
            Threads.@spawn begin
                rng = MersenneTwister(10*N+stepT)
                T = Temps[stepT]
                verbose && println("| Process strarted on thread #$(Threads.threadid()): T = $(T)")
                if from_infinity
                    spins = rand([1, -1], (N, N))
                else
                    spins = ones(Int64, (N, N))
                end
    
                P_add = isingwolff_Padd(T)
                e1 = @elapsed for i in 1:eqsteps
                    isingwolff_step!(N, spins, P_add; rng=rng)
                end
                τ = autocorr_times[stepT]
                e2 = @elapsed begin
                    uncorrelated_spins = ising_getuncorrconfigs!(N, spins, T, τ, n_uncorr; ntau=ntau, rng=rng) 
                end
    
                filename="ising_uncorr_configs_temp$(T)_size$(N).txt"
                open(joinpath([szpath, filename]), mode) do io
                    writedlm(io, reshape(uncorrelated_spins, (N*N, n_uncorr)), ',')
                end;
                verbose && println("| Process complete on thread #$(Threads.threadid()) (T = $T) in $(round(e1+e2, digits=3)) seconds.") 
            end
        end
        verbose && println("| Done")
        verbose && println(".==================================")
    end
end

function ising_getuncorrconfigs!(N, spins::Matrix, T, τ, n_uncorr; ntau=5, rng=TaskLocalRNG())
    ntau_τ = ntau*τ
    nsteps = ntau_τ*n_uncorr
    uncorrelated_spins = zeros(Float64, (N, N, n_uncorr))
    P_add = isingwolff_Padd(T)
    for j=1:nsteps
        isingwolff_step!(N, spins, P_add; rng=rng)
        if j%ntau_τ == 0
            uncorrelated_spins[:, :, j÷ntau_τ] = spins
        end
    end
    return uncorrelated_spins
end

function xy_getconfigdata_to_txt(
    lattice_sizes::AbstractArray{Int64},
    Temps::AbstractArray{Float64},
    n_uncorr::Int64,
    eqsteps::Int64; 
    autocorr_times::AbstractArray{Int64}=ones(Int64, length(Temps)),
    store_at::AbstractString="",
    verbose=true,
    mode="w",
    ntau=2
    )
    for N in lattice_sizes
        verbose && println(".==================================")
        verbose && println("| Lattice Size: $(N) x $(N)        ")
        verbose && println(".==================================")
        verbose && println("|  ")

        szpath = joinpath([store_at, "xy", "uncorr_configs", "Size$N"])
        ispath(szpath) ? 1 : mkpath(szpath)
        @sync for stepT in 1:length(Temps)
            Threads.@spawn begin
                rng = MersenneTwister(10*N+stepT)
                T = Temps[stepT]
                verbose && println("| Process strarted on thread #$(Threads.threadid()) (T = $(T)).")
        
                spins = zeros(Float64, (N, N))
                e1 = @elapsed for i in 1:eqsteps
                    xywolff_step!(N, spins, T; rng=rng)
                end

                τ = autocorr_times[stepT]
                e2 = @elapsed begin
                    uncorrelated_spins = xy_getuncorrconfigs!(N, spins, T, τ, n_uncorr; ntau=ntau, rng=rng)
                end
                
                filename="xy_uncorr_configs_temp$(T)_size$(N).txt"
                open(joinpath([szpath, filename]), mode) do io
                    writedlm(io, reshape(uncorrelated_spins, (N*N, n_uncorr)), ',')
                end;
                verbose && println("| Process complete on thread #$(Threads.threadid()) (T = $T) in $(round(e1+e2, digits=3)) seconds.")
            end
        end
        verbose && println("| Done.")
        verbose && println(".==================================")
    end
    nothing
end

function xy_getuncorrconfigs!(N, spins::Matrix, T, τ, n_uncorr; ntau=5, rng=TaskLocalRNG())
    ntau_τ = ntau*τ
    nsteps = ntau_τ*n_uncorr
    uncorrelated_spins = zeros(Float64, (N, N, n_uncorr))
    for j=1:nsteps
        xywolff_step!(N, spins, T; rng=rng)
        if j%ntau_τ == 0
            uncorrelated_spins[:, :, j÷ntau_τ] = spins
        end
    end
    return uncorrelated_spins
end

function ising_getcorrtime(
    N::Int64, T::Float64, msteps=6000, wsteps=250;
    wolff=false, from_infinity::Bool=false, eqsteps=1000
    )
    if from_infinity
        spins = rand([1.0, -1.0], (N, N))
    else
        spins = fill(1.0, (N, N))
    end
    mags = zeros(Float64, msteps)
    if wolff
        P_add = isingwolff_Padd(T)
        for i in 1:eqsteps
            isingwolff_step!(N, spins, P_add)
        end
        mags[1] = ising_total_magnetization(spins)
        for i in 1:msteps-1
            isingwolff_step!(N, spins, P_add)
            mags[i+1] = ising_total_magnetization(spins)
        end
    else
        E0 = ising_total_energy(spins)
        mags[1] = ising_total_magnetization(spins)
        for i in 1:eqsteps
            E0, mags[1] = isingmetro_step!(spins, T, E0, mags[i])
        end
        for i in 1:msteps-1
            E0, mags[i+1] = isingmetro_step!(spins, T, E0, mags[i])
        end
    end
    
    @. mags /= N^2
    corrfn = autocorrelation_fn(mags)

    # Measure the integrated correlation time in a small window
    τ = corrfn[1:wsteps] |> sum |> ceil
    return convert(Int64, τ)
end

function ising_getcorrtime(
    N::Int64, Temps::Vector{Float64}, msteps=50000, wsteps=1000;
    wolff=false, from_infinity::Bool=false, verbose::Bool=false, eqsteps=1000
    )
    corr_times = zeros(Int64, length(Temps))
    @sync for i=eachindex(Temps)
        Threads.@spawn begin
            T = Temps[i]
            verbose && println("> Process strarted on thread #$(Threads.threadid()): T = $(T)")
            corr_times[i] = ising_getcorrtime(N, T, msteps, wsteps; wolff=wolff, from_infinity=from_infinity, eqsteps=eqsteps)
            verbose && println("> Process complete on thread #$(Threads.threadid()): T = $T, τ=$(corr_times[i])")
        end
    end
    verbose && println("Results: $corr_times")
    return corr_times
end

function xy_getcorrtime!(spins::Matrix, T, msteps=20000)
    N = size(spins)[1]
    ergs = zeros(Float64, msteps)
    for i in 1:msteps-1
        xywolff_step!(spins, N, T)
        ergs[i+1] = xy_total_energy(spins, N)
    end
    
    @. ergs /= N^2
    corrfn = autocorrelation_fn(ergs)

    # Measure the integrated correlation time in a small window
    τ = corrfn[1:500] |> sum |> ceil
    return convert(Int64, τ)
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