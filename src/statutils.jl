"""
    autocorrelation_fn(mags, N)

Calculate the autocorrelation function (normalized) of the given time series array.
"""
function autocorrelation_fn(series)
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
function bootstrap_err(samples, calc_qty, args...; r = 200)
    nob = length(samples)
    resample_arr = zeros(Float64, r)
    for i in eachindex(resample_arr)
        resample = rand(samples, nob)
        resample_arr[i] = calc_qty(resample, args...)
    end
    err = std(resample_arr, corrected = false)
    return err
end


"""
    blocking_err(samples, calc_qty; blocks=20)

Estimate the error in the given samples by blocking method.
Here, `calc_qty` is the function to calculate the quantity in which error has to be calculated.
And, `blocks` is a keyword arguments giving number of blocks.
"""
function blocking_err(samples, calc_qty, args...; blocks = 20)
    block_array = zeros(Float64, blocks)
    blocklength = length(samples) ÷ blocks
    for i = 1:blocks
        sample_block = samples[(i-1)*blocklength+1:i*blocklength]
        block_array[i] = calc_qty(sample_block, args...)
    end
    err = std(block_array)
    return err
end

"""
Calculate site-site correlation function of given Matrix
"""
function ss_correlation_fn(sites::Matrix, N::Int64; metric=*)
    ss_corrs = zeros(Float64, N)
    nsamples = zeros(Float64, N)
    for i=1:N
        for j=1:N
            r = abs(i-j)
            ss_corrs[r+1] += metric(sites[i,i], sites[i,j]) 
            ss_corrs[r+1] += metric(sites[i,i], sites[j,i])
            nsamples[r+1] += 2
        end
    end
    return ss_corrs ./ nsamples
end

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
