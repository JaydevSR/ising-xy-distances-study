"""
    isingwolff_step!(spins, P_add)

Performs one cluster flip of the Ising lattice `spins` with selection probability `P_add`.
"""
function isingwolff_step!(N, spins, P_add)
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
            if spins[nn...] == sval && !cluster[nn...] && rand() < P_add
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
    isingwolff_Padd(T; J=1)

Calculate the probability of adding a neighbour to a cluster at temperature `T` and interaction energy `J`.
"""
function isingwolff_Padd(T; J=1)
    return 1 - exp(-2 * J / T)
end