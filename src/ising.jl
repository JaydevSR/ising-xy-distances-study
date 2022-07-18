# Model
mutable struct ClassicalIsing2D{T, LL}
    L::Int
    lattice::Matrix{T}
    counts::SVector{Int, 2}
    shifts::SVector
    # TODO
    # H::Float64
    # M::Float64
end

function ClassicalIsing2D(L::Int, start::Symbol=:cold)
    if start == :cold
        lattice = ones(Int, L, L)
        counts=(1, 0)
    elseif start == :hot
        lattice = rand([-1, 1], (L, L))
        count_1 = 0
        for i in eachindex(lattice)
            if lattice[i] == 1
                count_1 += 1
            end
        end
        counts = @SVector [count_1, L*L-count_1]
    else
        error("Start state can be one of symbols :$(:cold) or :$(:hot)")
    end

    # Square lattice with 4 nearest neighbors
    shifts = SA[
        CartesianIndex(1, 0), CartesianIndex(L - 1, 0), 
        CartesianIndex(0, 1), CartesianIndex(0, L - 1)
        ]
    return ClassicalIsing2D{Int, L}(L, lattice, counts, shifts)
end

function get_neighbors(model::ClassicalIsing2D, k::CartesianIndex)
    return SA[[CartesianIndex(mod1.((k+δ).I, model.L)) for δ in model.shifts]...]
end

magnetization(model::ClassicalIsing2D) = (model.counts[1] - model.counts[2])

function hamiltonian(model::ClassicalIsing2D)
    running_sum = 0
    for site in eachindex(model.lattice)
        for nn in get_neighbors(model, site)
            @inbounds running_sum += model.lattice[site] * model.lattice[nn]
        end
    end
    return - running_sum / 2  # divide by 2 because each bond counted twice
end

"""
    isingmetro_step!(model, T, E, M)

Perform one sweep of the lattice using single-spin-flip dynamics (1 sweep == N*N flip attempts).
Here arguments E and M are the total energy and total magnetization before the sweep.
Returns total energy and magnetization after sweep.
"""
function metropolis_update!(model, T, E, M)
    # TODO: change counts
    for i = 1:model.L^2
        k = rand(1:model.L, 2)
        ΔE = ising_delE_flip(model, k)
        if isingmetro_accept(ΔE, T)
            model.lattice[k] *= -1
            E = E + ΔE
            M = M + 2model.lattice[k]
        end
    end
    return E, M
end


"""
    isingmetro_accept(spins, N, k_i, k_j, T)

Determine whether to accept the next state or not according to Metropolis acceptance rates.
Returns `true` or `false`.
"""
@inline function isingmetro_accept(ΔE, T)
    # Metropolis Acceptance Rates
    if ΔE <= 0
        return true
    elseif rand() < exp(-ΔE / T)
        return true
    else
        return false
    end
end


"""
    ising_delE_flip(model, k)

Calculate the energy difference between two states for one spin flip at site `k`.
"""
@inline function ising_delE_flip(model, k)
    ΔE = 0
    for nn in get_neighbors(model, k)
        ΔE += model.lattice[nn]
    end
    ΔE *= 2model.lattice[k]
end

"""
    wolff_cluster_update(model, P_add; [rng=TaskLocalRNG(), stack=LazyStack()])

Performs one cluster flip of the Ising lattice `spins` at temperature `T`.
"""
function wolff_update!(
                    model::ClassicalIsing2D,
                    P_add::Float64;
                    rng=TaskLocalRNG(),
                    stack::LazyStack=LazyStack(CartesianIndex{2}))
    cluster = falses(model.L, model.L)
    seed = CartesianIndex(Tuple(rand(1:model.L, 2)))
    empty!(stack)
    push!(stack, seed)
    @inbounds sval = model.lattice[seed]
    @inbounds cluster[seed] = true
    n_flips = 0
    while !isempty(stack)
        k = pop!(stack)
        @inbounds model.lattice[k] *= -1
        n_flips += 1
        for nn in get_neighbors(model, k)
            if model.lattice[nn] == sval && !cluster[nn] && rand(rng) < P_add
                push!(stack, nn)
                @inbounds cluster[nn] = true
            end
        end
    end
    if sval == 1
        model.counts = @SVector [model.counts[1] - n_flips, model.counts[2] + n_flips]
    else
        model.counts = @SVector [model.counts[1] + n_flips, model.counts[2] - n_flips]
    end
    nothing
end

"""
Probability of adding a site to cluster
"""
isingwolff_Padd(T::Float64) = 1 - exp(-2 / T)