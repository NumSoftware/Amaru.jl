mutable struct Mat6x6<:AbstractArray{Float64,2}
    data::NTuple{36, Float64}

    function Mat6x6()
        return new()
    end

    function Mat6x6(vals::Real...)
        return new((vals...,))
    end

    function Mat6x6(data::NTuple{36, Float64})
        return new(data)
    end
    function Mat6x6(M::Mat6x6)
        return new(M.data)
    end
    function Mat6x6(M::Array{<:Real,2})
        return new(tuple(M...))
    end
end

Base.convert(::Type{Mat6x6}, A::AbstractArray{Float64,2}) = Mat6x6(A)

macro Mat6x6(expr)
    @assert expr.head==:vcat
    @assert length(expr.args)==6
    @assert length(expr.args[1].args)==6
    args = [ expr.args[i].args[j] for j in 1:6 for i in 1:6 ]
    return esc(Expr(:call, :Mat6x6, args...))
end


Base.length(::Mat6x6) = 36
Base.size(::Mat6x6) = (6,6)
Base.strides(::Mat6x6) = (1,6)
Base.unsafe_convert(::Type{Ptr{Float64}}, A::Mat6x6) = Base.unsafe_convert(Ptr{Float64}, pointer_from_objref(A))
Base.zeros(::Type{Mat6x6}) = Mat6x6(0.0.*Mat6x6().data)
Base.copy(M::Mat6x6) = Mat6x6(M.data)


function Base.getindex(M::Mat6x6, k::Int)
    0<k<=36 || throw(BoundsError(M,i))
    return M.data[k]
    # return GC.@preserve M unsafe_load(Base.unsafe_convert(Ptr{Float64}, pointer_from_objref(M)), k)
end

function Base.getindex(M::Mat6x6, i::Int, j::Int)
    0<i<=6 || throw(BoundsError(M,(i,j)))
    0<j<=6 || throw(BoundsError(M,(i,j)))
    
    k = 6*(j-1) + i
    return M.data[k]
end


function Base.setindex!(M::Mat6x6, val::Float64, i::Int, j::Int)
    0<i<=6 || throw(BoundsError(M,(i,j)))
    0<j<=6 || throw(BoundsError(M,(i,j)))

    k = 6*(j-1) + i
    GC.@preserve M unsafe_store!(Base.unsafe_convert(Ptr{Float64}, pointer_from_objref(M)), val, k)
end


Base.:*(a::Real, M::Mat6x6) = Mat6x6(a.*M.data)
Base.:*(M::Mat6x6, a::Real) = Mat6x6(a.*M.data)

function Base.:*(M::Mat6x6, X::Vec6)
    d = M.data
    Vec6( 
        d[1]*X.v1 + d[07]*X.v2 + d[13]*X.v3 + d[19]*X.v4 + d[25]*X.v5 + d[31]*X.v6,
        d[2]*X.v1 + d[08]*X.v2 + d[14]*X.v3 + d[20]*X.v4 + d[26]*X.v5 + d[32]*X.v6,
        d[3]*X.v1 + d[09]*X.v2 + d[15]*X.v3 + d[21]*X.v4 + d[27]*X.v5 + d[33]*X.v6,
        d[4]*X.v1 + d[10]*X.v2 + d[16]*X.v3 + d[22]*X.v4 + d[28]*X.v5 + d[34]*X.v6,
        d[5]*X.v1 + d[11]*X.v2 + d[17]*X.v3 + d[23]*X.v4 + d[29]*X.v5 + d[35]*X.v6,
        d[6]*X.v1 + d[12]*X.v2 + d[18]*X.v3 + d[24]*X.v4 + d[30]*X.v5 + d[36]*X.v6,
    )
end

Base.:*(M::Mat6x6, X::AbstractArray{Float64,1}) = return M*Vec6(X)
Base.:+(M::Mat6x6, N::Mat6x6) = return Mat6x6(M.data.+N.data)
Base.:-(M::Mat6x6, N::Mat6x6) = return Mat6x6(M.data.+N.data)
Base.:+(M::Mat6x6, N::AbstractArray{Float64,2}) = return M + Mat6x6(N)
Base.:+(M::AbstractArray{Float64,2}, N::Mat6x6) = return Mat6x6(M) + N