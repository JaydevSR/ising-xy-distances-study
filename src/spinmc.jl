using Statistics
using DelimitedFiles
using Random
using StaticArrays

include("lazystack.jl")
include("ising.jl")
inclure("statutils.jl")

const ising_Tc = 2 / log1p(sqrt(2))
