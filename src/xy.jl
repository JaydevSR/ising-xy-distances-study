# Model
mutable struct ClassicalXYModel2D{T, LL} <: SpinModel2D{T}
    L::Int
    lattice::Matrix{T}
    proj_X::Float64
    proj_Y::Float64
    shifts::SVector
end

function hamiltonian(model::ClassicalXYModel2D)
    running_sum = 0
    for idx in CartesianIndices(model.lattice)
        s_k = model.lattice[idx]
        for nn in get_neighbors(model, idx)
            running_sum += cos2pi(s_k - model.lattice[nn])
        end
    end
    return -running_sum / 2  # divide by 2 because each bond counted twice
end

function structure_factor(model::ClassicalXYModel2D; scaled::Bool=false)
    R_spins = 0
    for i in CartesianIndices(model.lattice)
        for j in CartesianIndices(model.lattice)
            R_spins += cos2pi(model.lattice[i] - model.lattice[j])
        end
    end
    if scaled
        R_spins /= N^2
    else
        R_spins /= N^4
    end
    return R_spins
end

# TODO: Update projections in the model
function wolff_update!(
                    model::ClassicalXYModel2D,
                    T::Float64;
                    rng=TaskLocalRNG(),
                    stack::LazyStack=LazyStack(CartesianIndex{2}))
    cluster = falses(model.L, model.L)
    seed = CartesianIndex(Tuple(rand(1:model.L, 2)))  # seed spin position
    u_flip = rand(rng)  # Random unit vector in xy plane
    push!(stack, seed)
    @inbounds cluster[seed] = true
    while !isempty(stack)
        k = pop!(stack)
        @inbounds kval = model.lattice[k]
        xywolff_flip_spin!(model, k, u_flip)
        for nn in get_neighbors(model, k)
            nnval = spins[nn]
            if !cluster[nn] && rand(rng) < xywolff_Padd(u_flip, nnval, kval, T)
                push!(stack, nn)
                @inbounds cluster[nn] = true
            end
        end
    end
end

@inline @inbounds function xywolff_flip_spin!(model, pos::CartesianIndex, u_flip::Float64)
    old = model.lattice[pos]
    new = 0.5 + 2 * u_flip - old  # flipping w.r.t vector with angle ϕ: θ --> π + 2ϕ - θ
    new = mod(new + 1, 1)
    model.lattice[pos] = new
    return model
end

function xywolff_Padd(u_flip::Float64, s1::Float64, s2::Float64, T::Float64)
    arg = -2*cos2pi(u_flip - s1)*cos2pi(u_flip - s2) / T
    return 1 - exp(arg)
end