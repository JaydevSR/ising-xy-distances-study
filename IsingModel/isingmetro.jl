# Functions

"""
    ising_metropolis_sweep(spins, T, E, M)

Perform one sweep of the lattice using single-spin-flip dynamics (`1 sweep == N*N flip attempts`).
Here arguments `E` and `M` are the total energy and total magnetization before the sweep.
Returns total energy and magnetization after sweep.
"""
function ising_metropolis_sweep!(spins, T, E, M)
    N = size(spins)[1]
    for i=1:N^2
        k = rand(1:N, 2)
        ΔE = deltaE(spins, k, N)
        if accept_flip(ΔE, T)
            spins[k[1], k[2]] *= -1
            E = E + ΔE
            M = M + 2spins[k[1], k[2]]
        end
    end
    return E, M
end

"""
    accept_flip(spins, N, k_i, k_j, T)

Determine whether to accept the next state or not according to Metropolis acceptance rates.
Returns `true` or `false`.
"""
function accept_flip(ΔE, T)
    # Metropolis Acceptance Rates
    if ΔE <= 0
        return true
    elseif rand() < exp(-ΔE/T)
        return true
    else
        return false
    end
end

"""
    delta_E(spins, k)

Calculate the energy difference between two states for one spin flip at site `k`.
"""
function deltaE(spins, k, N)
    ΔE = 0
    for δ ∈ ([1, 0], [N-1, 0], [0, 1], [0, N-1])
        nn = k + δ
        @. nn = mod1(nn, N)  # Apply periodic boundary conditions
        ΔE += spins[nn[1], nn[2]]
    end
    ΔE *= 2spins[k[1], k[2]]
end

"""
    total_magnetization(spins)

Calculate the total magnetization of square spin lattice.
"""
function total_magnetization(spins)
    return sum(sum(spins))
end

"""
    total_energy(spins)

Calculate the total energy of the square spin lattice (with zero field and J=1).
"""
function total_energy(spins)
    N = size(spins)[1]
    running_sum = 0
    for i=1:N
        for j=1:N
            s_k = spins[i, j]
            for δ ∈ ([1, 0], [N-1, 0], [0, 1], [0, N-1])
                nn = [i, j] + δ
                @. nn = mod1(nn, N)  # Apply periodic boundary conditions
                running_sum += s_k*spins[nn[1], nn[2]]
            end
        end
    end
    return -running_sum / 2  # divide by 2 because each bond counted twice
end

"""
    autocorrelation_fn(mags, N)

Calculate the autocorrelation function (normalized) of the given time series array.
"""
function autocorrelation_fn(series, N)
    tmax = length(series)
    autocorr = zeros(Float64, tmax)
    for t ∈ 1:tmax-1
        sum1 = 0
        sum2 = 0
        sum3 = 0
        for tk ∈ 1:tmax-t
            sum1 += series[tk]*series[tk+t]
            sum2 += series[tk]
            sum3 += series[tk+t]
        end
        autocorr[t] = sum1 / (tmax-t) - (sum2*sum3) / (tmax - t)^2
    end
    @. autocorr /= N^2
    @. autocorr /= autocorr[1]
    return autocorr
end