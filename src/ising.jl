# Model
mutable struct ClassicalIsingModel2D{T, LL} <: SpinModel2D{T}
    L::Int
    lattice::Matrix{T}
    counts::SVector{2, Int}
    shifts::SVector
end

const ising_Tc = 2 / log1p(sqrt(2))

function ClassicalIsingModel2D(L::Int, start::Symbol=:cold)
    if start == :cold
        lattice = ones(Int, L, L)
        counts=(L*L, 0)
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
    return ClassicalIsingModel2D{Int, L}(L, lattice, counts, shifts)
end

@inbounds function magnetization(model::ClassicalIsingModel2D)
    m = (model.counts[1] - model.counts[2])
    return abs(m)
end

function hamiltonian(model::ClassicalIsingModel2D)
    running_sum = 0
    for site in eachindex(model.lattice)
        for nn in get_neighbors(model, site)
            @inbounds running_sum += model.lattice[site] * model.lattice[nn]
        end
    end
    return - running_sum / 2  # divide by 2 because each bond counted twice
end

# TODO: simplify this using counts
function structure_factor(model::ClassicalIsingModel2D; scaled::Bool=false)
    R_spins = 0
    for i in CartesianIndices(model.lattice)
        for j in CartesianIndices(model.lattice)
            R_spins += xy_spindot(model.lattice[i], model.lattice[j])
        end
    end
    if scaled
        R_spins /= N^2
    else
        R_spins /= N^4
    end
    return R_spins
end

# Algorithms

function wolff_update!(
                    model::ClassicalIsingModel2D,
                    T::Float64,
                    cluster::BitMatrix=falses(model.L, model.L),
                    stack::LazyStack=LazyStack(CartesianIndex{2});
                    P_add::Float64=isingwolff_Padd(T),
                    )
    empty!(stack)
    cluster .= false
    seed = CartesianIndex(Tuple(rand(1:model.L, 2)))
    push!(stack, seed)
    @inbounds sval = model.lattice[seed]
    @inbounds cluster[seed] = true
    n_flips = 0
    nn = CartesianIndex(0, 0)
    while !isempty(stack)
        k = pop!(stack)
        @inbounds model.lattice[k] *= -1
        n_flips += 1
        for δ in model.shifts
            @inbounds nn = CartesianIndex(mod1(k[1] + δ[1], model.L), mod1(k[2] + δ[2], model.L))
            if model.lattice[nn] == sval && !cluster[nn] && rand() < P_add
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
    return model
end

isingwolff_Padd(T::Float64) = 1 - exp(-2 / T)