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
        Threads.@threads for stepT in 1:length(Temps)
            T = Temps[stepT]
            verbose && println("| Process strarted on thread #$(Threads.threadid()): T = $(T)")
            if from_infinity
                spins = rand([1.0, -1.0], (N, N))
            else
                spins = fill(1.0, (N, N))
            end

            P_add = isingwolff_Padd(T)
            e1 = @elapsed for i in 1:eqsteps
                isingwolff_step!(N, spins, P_add)
            end
            τ = autocorr_times[stepT]
            e2 = @elapsed begin
                uncorrelated_spins = ising_getuncorrconfigs!(N, spins, T, τ, n_uncorr; ntau=ntau) 
            end

            filename="ising_uncorr_configs_temp$(T)_size$(N).txt"
            open(joinpath([szpath, filename]), mode) do io
                writedlm(io, reshape(uncorrelated_spins, (N*N, n_uncorr)), ',')
            end;
            verbose && println("| Process complete on thread #$(Threads.threadid()) (T = $T) in $(round(e1+e2, digits=3)) seconds.")
        end
        verbose && println("| Done.")
        verbose && println(".==================================")
    end
end

function ising_getuncorrconfigs!(N, spins::Matrix, T, τ, n_uncorr; ntau=5)
    ntau_τ = ntau*τ
    nsteps = ntau_τ*n_uncorr
    uncorrelated_spins = zeros(Float64, (N, N, n_uncorr))
    P_add = isingwolff_Padd(T)
    for j=1:nsteps
        isingwolff_step!(N, spins, P_add)
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
        for stepT in 1:length(Temps)
            T = Temps[stepT]
            verbose && println("| Process strarted on thread #$(Threads.threadid()) (T = $(T)).")
    
            spins = zeros(Float64, (N, N))
            e1 = @elapsed for i in 1:eqsteps
                xywolff_step!(N, spins, T)
            end

            τ = autocorr_times[stepT]
            e2 = @elapsed begin
                uncorrelated_spins = xy_getuncorrconfigs!(N, spins, T, τ, n_uncorr; ntau=ntau)
            end
            
            filename="xy_uncorr_configs_temp$(T)_size$(N).txt"
            open(joinpath([szpath, filename]), mode) do io
                writedlm(io, reshape(uncorrelated_spins, (N*N, n_uncorr)), ',')
            end;
            verbose && println("| Process complete on thread #$(Threads.threadid()) (T = $T) in $(round(e1+e2, digits=3)) seconds.")
        end
        verbose && println("| Done.")
        verbose && println(".==================================")
    end
    nothing
end

function xy_getuncorrconfigs!(N, spins::Matrix, T, τ, n_uncorr; ntau=5)
    ntau_τ = ntau*τ
    nsteps = ntau_τ*n_uncorr
    uncorrelated_spins = zeros(Float64, (N, N, n_uncorr))
    for j=1:nsteps
        xywolff_step!(N, spins, T)
        if j%ntau_τ == 0
            uncorrelated_spins[:, :, j÷ntau_τ] = spins
        end
    end
    return uncorrelated_spins
end
