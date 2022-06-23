"""
Calculate structure factor of given Matrix of spins
"""
function structure_factor(spins::Matrix, N::Int64; metric=*, scaled::Bool=false)
    R_spins = 0
    for i in CartesianIndices(spins)
        for j in CartesianIndices(spins)
            R_spins += metric(spins[i], spins[j])
        end
    end
    if scaled
        R_spins /= N^2
    else
        R_spins /= N^4
    end
    return R_spins
end
