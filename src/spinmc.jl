using Statistics
using DelimitedFiles
using Random
using StaticArrays

include("lazystack.jl")

abstract type SpinModel2D{T} end

function get_neighbors(model::SpinModel2D, k::CartesianIndex)
    ns = length(model.shifts)
    return ntuple(i -> CartesianIndex(mod1.((k+model.shifts[i]).I, model.L)), ns)
end

include("ising.jl")

include("xy.jl")

function cos2pi(x::Float64)
    return cospi(2x)
end

function sin2pi(x::Float64)
    return sinpi(2x)
end

inclure("statutils.jl")
