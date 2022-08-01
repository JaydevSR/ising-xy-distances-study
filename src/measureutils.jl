function ising_get_measurements_to_txt(
    lattice_sizes::AbstractArray{Int64},
    temps::AbstractArray{Float64},
    nconfigs::Int64,
    eqsteps::Int64,
    autocorr_times::AbstractArray{Int64}=ones(Int64, length(Temps));
    ntau=2,
    store_at::AbstractString="",
    start::Symbol=:cold,
    verbose=true,
    mode="w",
    get_mags::Bool=true,
    get_structure_factors::Bool=true
    )
    for L in lattice_sizes
        verbose && println(".==================================")
        verbose && println("| Lattice Size: $(L) x $(L)        ")
        verbose && println(".==================================")
        verbose && println("|  ")
        @sync for stepT in 1:length(temps)
            Threads.@spawn ising_get_measurements_to_txt(
                L, $temps[stepT], nconfigs, eqsteps, $autocorr_times[stepT];
                store_at=store_at, start=start, ntau=ntau,
                mode=mode, verbose=verbose, get_mags=get_mags, 
                get_structure_factors=get_structure_factors
                )
        end
        verbose && println("| Done.")
        verbose && println(".==================================")
    end
    nothing
end

function ising_get_measurements_to_txt(
    L::Int,
    T::Float64,
    nconfigs::Int64,
    eqsteps::Int64,
    τ::Int64=1;
    ntau=2,
    start::Symbol=:cold,
    store_at::AbstractString="",
    verbose=true,
    mode="w",
    get_mags::Bool=true,
    get_structure_factors::Bool=true
    )
    verbose && println("| Process strarted for T = $(T).")
    el = @elapsed begin
        model = ClassicalIsingModel2D(L, start)
        szpath = joinpath([store_at, "Ising", "Size$L"])
        mkpath(szpath)
        P_add = isingwolff_Padd(T)
        stack = LazyStack(Int)
        cluster = falses(L, L)
        el_eq = @elapsed for i=1:eqsteps  # equilibration
            wolff_update!(model, T; P_add=P_add, stack=stack, cluster=cluster)
        end
    
        uncorrelated_spins = zeros(Int64, (model.L^2, nconfigs))
        get_mags && (mags = zeros(Float64, nconfigs))
        get_structure_factors && (struc_facs = zeros(Float64, nconfigs))

        nτ = ntau*τ
        numsteps = nτ*nconfigs
        el = @elapsed for j in 1:numsteps
            wolff_update!(model, T; P_add=P_add, stack=stack, cluster=cluster)
            if j%nτ == 0
                uncorrelated_spins[:, j÷nτ] = reshape(model.lattice, model.L^2)
                get_mags && (mags[j÷nτ] = magnetization(model))
                get_structure_factors && (struc_facs[j÷nτ] = structure_factor(model))
            end
        end
    
        filename_conf="ising_uncorr_configs_temp$(T)_L$(L).txt"
        mkpath(joinpath([szpath, "uncorr_configs"]))
        open(joinpath([szpath, "uncorr_configs", filename_conf]), mode) do io
            writedlm(io, uncorrelated_spins, ',')
        end;
        if get_mags
            filename_mag="ising_mags_temp$(T)_L$(L).txt"
            mkpath(joinpath([szpath, "mags"]))
            open(joinpath([szpath, "mags", filename_mag]), mode) do io
                writedlm(io, mags, ',')
            end;
        end
        if get_structure_factors
            filename_struc_fac="ising_struc_facs_temp$(T)_L$(L).txt"
            mkpath(joinpath([szpath, "struc_facs"]))
            open(joinpath([szpath, "struc_facs", filename_struc_fac]), mode) do io
                writedlm(io, struc_facs, ',')
            end;
        end
    end
    verbose && println("| Process complete for T = $T in $(round(el+el_eq, digits=0)) seconds.")
    nothing
end

function ising_getcorrtime(
    L::Int64, T::Float64, msteps=20000, wsteps=1000, eqsteps=1000; verbose=true)
    verbose && (@info "Calculating autocorrelation time for system size $(L)x$(L) and temperature $(T).")
    model = ClassicalIsingModel2D(L, :cold)
    mags = zeros(Float64, msteps)
    P_add = isingwolff_Padd(T)
    stack = LazyStack(CartesianIndex{2})
    cluster = falses(L, L)
    verbose && (@info "Equilibrating the system over $(eqsteps) steps.")
    for i=1:eqsteps  # equilibration
        wolff_update!(model, T; P_add=P_add, stack=stack, cluster=cluster)
    end

    verbose && (@info "Calculating magnetizations over $(msteps) steps.")
    for i=1:msteps
        wolff_update!(model, T; P_add=P_add, stack=stack, cluster=cluster)
        @inbounds mags[i] = magnetization(model)
    end
    
    verbose && (@info "Calculating autocorrelation function.")
    @. mags /= L^2
    corrfn = autocorrelation_fn(mags)

    # Measure the integrated correlation time in a small window
    @views τ = corrfn[1:wsteps] |> sum |> ceil
    τ = convert(Int64, τ)
    verbose && (@info "Autocorrelation time: $(τ)")
    return τ
end