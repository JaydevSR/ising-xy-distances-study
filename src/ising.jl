"""
    isingmetro_step!(spins, T, E, M)

Perform one sweep of the lattice using single-spin-flip dynamics (1 sweep == N*N flip attempts).
Here arguments E and M are the total energy and total magnetization before the sweep.
Returns total energy and magnetization after sweep.
"""
function isingmetro_step!(spins, T, E, M)
    N = size(spins)[1]
    for i = 1:N^2
        k = rand(1:N, 2)
        ΔE = ising_delE_flip(spins, k, N)
        if isingmetro_accept(ΔE, T)
            spins[k[1], k[2]] *= -1
            E = E + ΔE
            M = M + 2spins[k[1], k[2]]
        end
    end
    return E, M
end


"""
    isingmetro_accept(spins, N, k_i, k_j, T)

Determine whether to accept the next state or not according to Metropolis acceptance rates.
Returns `true` or `false`.
"""
function isingmetro_accept(ΔE, T)
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
    ising_delE_flip(spins, k)

Calculate the energy difference between two states for one spin flip at site k.
"""
function ising_delE_flip(spins, k, N)
    ΔE = 0
    for δ ∈ ([1, 0], [N - 1, 0], [0, 1], [0, N - 1])
        nn = k + δ
        @. nn = mod1(nn, N)  # Apply periodic boundary conditions
        ΔE += spins[nn[1], nn[2]]
    end
    ΔE *= 2spins[k[1], k[2]]
end

"""
    isingwolff_step!(spins, P_add)

Performs one cluster flip of the Ising lattice `spins` with selection probability `P_add`.
"""
function isingwolff_step!(N, spins, P_add; rng=TaskLocalRNG())
    nbrs = [[1, 0], [N - 1, 0],
            [0, 1], [0, N - 1]]
    seed = rand(1:N, 2)  # seed spin position
    stack = []
    sizehint!(stack, N*N)
    push!(stack, seed)
    cluster = falses(size(spins))
    @inbounds sval = spins[seed...]
    @inbounds cluster[seed...] = true
    while !isempty(stack)
        k = pop!(stack)
        @inbounds spins[k...] = -spins[k...]
        for δ ∈ nbrs
            nn = k + δ
            @. nn = mod1(nn, N)  # Apply periodic boundary conditions
            if spins[nn...] == sval && !cluster[nn...] && rand(rng) < P_add
                push!(stack, nn)
                @inbounds cluster[nn...] = true
            end
        end
    end
    # ΔM = 2 * spins[seed...] * sum(cluster)
    # return ΔM
    nothing
end

"""
    isingwolff_get_cluster!(N, spins, seed, P_add)

Generates a cluster in the Ising lattice `spins` at site `seed` with neighbour selection probability `P_add`.
"""
function isingwolff_get_cluster!(N, spins, seed, P_add, nbrs)
    cluster = falses(size(spins))
    @inbounds sval = spins[seed...]
    stack = [seed]
    @inbounds cluster[seed...] = true
    while !isempty(stack)
        k = pop!(stack)
        for δ ∈ nbrs
            nn = k + δ
            @. nn = mod1(nn, N)  # Apply periodic boundary conditions
            if spins[nn...] == sval && !cluster[nn...] && rand() < P_add
                push!(stack, nn)
                @inbounds cluster[nn...] = true
            end
        end
    end
    return cluster
end

"""
    isingwolff_Padd(T; J=1)

Calculate the probability of adding a neighbour to a cluster at temperature `T` and interaction energy `J`.
"""
function isingwolff_Padd(T; J=1)
    return 1 - exp(-2 * J / T)
end