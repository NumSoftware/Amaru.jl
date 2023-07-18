

# struct SArray{S, T, N, L} <: AbstractArray{T, N} where {S,L}
struct SArray{S, T, N, L} <: AbstractArray{T, N} where {S,L}
    data::NTuple{L,T}

    function SArray{S, T, N, L}(X::NTuple{L,T}) where {S, T, N, L}
        return new{S, T, N, L}(X)
    end
    
    function SArray{S, T, N, L}(X::AbstractArray) where {S, T, N, L}
        return new{size(X), eltype(X), length(X), prod(S)}( (X...,) )
    end
end

function SArray(A::Array{Float64,N}) where N
    S = size(A)
    L = length(A)
    return SArray{S,Float64,N,L}(tuple(A...))
end

function SArray{S, N, L}(X::Tuple) where {S,N,L}
    if all(isbits,X)
        Y = promote(X...)
        return SArray{S, eltype(Y), N, L}(Y)
    end
    return SArray{S, eltype(X), N, L}(X)
end

# Display function for 1D SArray
function Base.display(A::SArray{S, T, 1, L}) where {S, T, L}
    n = length(A)
    println("$n-element $(typeof(A)):")
    for i in 1:n
        println(" ", A[i])
    end
end

macro SArray(expr)
    @assert expr.head in (:vect, :vcat, :hcat)

    if expr.head in (:vect, :hcat)
        L = length(expr.args)
        S = expr.head==:vect ? (L,) : (1,L)
        args = ((e for e in expr.args)...,)
    elseif expr.head==:vcat
        n = length(expr.args) # rows
        m = length(expr.args[1].args) # cols
        S = (n, m)
        L = m*n
        args = ((expr.args[i].args[j] for j in 1:m for i in 1:n)...,)
    end

    if all(isbits, args)
        args = promote(args...)
    end

    N = length(S)
    T = eltype(args)

    if args isa NTuple && isbitstype(T)
        return quote
            SArray{$S, $T, $N, $L}(($(args...),))
        end
    else
        escargs = tuple( (esc(x) for x in args)..., )
        return quote
            SArray{$S, $N, $L}(($(escargs...),))
        end
    end
end

Base.length(::SArray{S,T,N,L}) where {S,T,N,L} = L
Base.size(::SArray{S,T,N,L}) where {S,T,N,L} = S

# Base.zeros(::Type{SArray}, dims::Int...) = SArray{dims, Float64, length(dims), prod(dims)}(ntuple(i->0.0, prod(dims)))

Base.strides(::SArray{S,T,N,L}) where {S,T,N,L} = Base.size_to_strides(1, S...)
Base.unsafe_convert(::Type{Ptr{T}}, A::SArray) where T = Base.unsafe_convert(Ptr{T}, pointer_from_objref(A))

function Base.getindex(A::SArray{S,T,N,L}, idxs::Int...) where {S,T,N,L}
    if length(idxs)==1
        k = idxs[1]
    elseif length(idxs)==2
        i, j = idxs
        k = S[1]*(j-1) + i
    end
    # @show idxs
    # @show k
    @assert k<=length(A.data)
    return A.data[k]
end