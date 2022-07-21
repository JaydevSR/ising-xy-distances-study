"""
    autocorrelation_fn(mags, N)

Calculate the autocorrelation function (normalized) of the given time series array.
"""
function autocorrelation_fn(series::Vector{<:Real})
    tmax = length(series)
    autocorr = zeros(Float64, tmax)
    for t ∈ 1:tmax-1
        sum1 = 0
        sum2 = 0
        sum3 = 0
        for tk ∈ 1:tmax-t
            sum1 += series[tk] * series[tk+t]
            sum2 += series[tk]
            sum3 += series[tk+t]
        end
        autocorr[t] = sum1 / (tmax - t) - (sum2 * sum3) / (tmax - t)^2
    end
    @. autocorr /= autocorr[1]
    return autocorr
end


"""
    bootstrap_err(samples, calc_qty; r=100)

Estimate the error in the given samples by bootstrap method.
Here, `calc_qty` is the function to calculate the quantity in which error has to be calculated.
And, `r` is a keyword arguments giving number of resamples.
"""
function bootstrap_err(samples::Vector{<:Real}, calc_qty::F, args...; r = 200) where F
    nob = length(samples)
    resample_arr = zeros(Float64, r)
    for i in eachindex(resample_arr)
        resample = rand(samples, nob)
        resample_arr[i] = calc_qty(resample, args...)
    end
    err = std(resample_arr, corrected = false)
    return err
end

function structure_factor(spins::Matrix, L::Int; metric=*, scaled::Bool=false)
    R_spins = 0
    for i in CartesianIndices(spins)
        for j in CartesianIndices(spins)
            R_spins += metric(spins[i], spins[j])
        end
    end
    if scaled
        R_spins /= L^2
    else
        R_spins /= L^4
    end
    return R_spins
end
