# Model
mutable struct ClassicalIsingModel2D{T, LL} <: SpinModel2D{T}
    L::Int
    lattice::Matrix{T}
    counts::NTuple{2, Int}
    shifts::NTuple{4, CartesianIndex{2}}
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
        counts = (count_1, L*L-count_1)
    else
        error("Start state can be one of symbols :$(:cold) or :$(:hot)")
    end

    # Square lattice with 4 nearest neighbors
    shifts = (
        CartesianIndex(1, 0), CartesianIndex(L - 1, 0), 
        CartesianIndex(0, 1), CartesianIndex(0, L - 1)
    )
    return ClassicalIsingModel2D{Int, L}(L, lattice, counts, shifts)
end

@inbounds function magnetization(model::ClassicalIsingModel2D)
    m = (model.counts[1] - model.counts[2])
    return abs(m)
end

function hamiltonian(model::ClassicalIsingModel2D)
    H = 0
    for site in eachindex(model.lattice)
        for nn in get_neighbors(model, site)
            @inbounds H = muladd(model.lattice[site], model.lattice[nn], H)
        end
    end
    return - H / 2  # divide by 2 because each bond counted twice
end

# TODO: simplify this using counts
function structure_factor(model::ClassicalIsingModel2D)
    R_spins = model.counts[1]*(model.counts[1] - model.counts[2])
    R_spins += model.counts[2]*(model.counts[2] - model.counts[1])
    R_spins /= model.L^4
    return R_spins
end

# Algorithms

function wolff_update!(
                    model::ClassicalIsingModel2D,
                    T::Float64;
                    cluster::BitMatrix=falses(model.L, model.L),
                    stack::LazyStack=LazyStack(Int),
                    P_add::Float64=isingwolff_Padd(T),
                    )
    if eltype(stack) != Int
        error("Stack must be of type Int")
    end
    empty!(stack)

    if size(cluster) != (model.L, model.L)
        error("Cluster must be of size $(model.L) x $(model.L)")
    end
    cluster .= false
    
    seedx, seedy = rand(1:model.L), rand(1:model.L)
    push!(stack, seedx)
    push!(stack, seedy)
    sval = model.lattice[seedx, seedy]
    cluster[seedx, seedy] = true
    n_flips = 0
    while !isempty(stack)
        ky = pop!(stack)  #! order matters
        kx = pop!(stack)
        model.lattice[kx, ky] *= -1
        n_flips += 1
        for δ in model.shifts
            nnx, nny = mod1(kx + δ[1], model.L), mod1(ky + δ[2], model.L)
            if model.lattice[nnx, nny] == sval && !cluster[nnx, nny] && rand() < P_add
                push!(stack, nnx)
                push!(stack, nny)
                cluster[nnx, nny] = true
            end
        end
    end
    if sval == 1
        model.counts = (model.counts[1] - n_flips, model.counts[2] + n_flips)
    else
        model.counts = (model.counts[1] + n_flips, model.counts[2] - n_flips)
    end
    return model
end

isingwolff_Padd(T::Float64) = -expm1(-2 / T) # 1 - exp(-2/T)