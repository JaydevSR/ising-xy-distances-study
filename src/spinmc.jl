using Statistics
using DelimitedFiles
using Random
using StaticArrays

include("lazystack.jl")

abstract type SpinModel2D{T} end

function get_neighbors(model::SpinModel2D, k::CartesianIndex)
    ns = length(model.shifts)
    return ntuple(i -> CartesianIndex(mod1(k[1] + model.shifts[i][1], model.L), mod1(k[2]+model.shifts[i][2], model.L)), ns)
end

include("ising.jl")

include("xy.jl")

function cos2pi(x::Float64)
    return cospi(2x)
end

function sin2pi(x::Float64)
    return sinpi(2x)
end

include("statutils.jl")
