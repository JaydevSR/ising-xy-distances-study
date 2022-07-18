mutable struct LazyStack{T}
    n::Int
    vals::Vector{T}
end

LazyStack() = LazyStack{Any}(0, [])
LazyStack(values::Vector{T}) where T <: Any = LazyStack{T}(length(values), values)
LazyStack(T::Type) = LazyStack{T}(0, [])

function Base.empty!(ls::LazyStack)
    ls.n = 0
    return ls
end

Base.isempty(ls::LazyStack) = (ls.n == 0)

Base.eltype(::LazyStack{T}) where T = T

function Base.push!(ls::LazyStack{T}, element::T) where T <: Any
    ls.n += 1
    if ls.n > length(ls.vals)
        push!(ls.vals, element)
    else
        ls.vals[ls.n] = element
    end
    return ls
end

function Base.pop!(ls::LazyStack{T}) where T <: Any
    if ls.n > 0
        ls.n -= 1
        return ls.vals[ls.n + 1]
    else
        error("Stack Underflow")
    end
end
