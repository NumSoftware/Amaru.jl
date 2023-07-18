mutable struct Vec6<:AbstractArray{Float64,1}
    v1::Float64
    v2::Float64
    v3::Float64
    v4::Float64
    v5::Float64
    v6::Float64

    function Vec6()
        return new()
    end
    function Vec6(v1, v2, v3, v4, v5, v6)
        return new(v1, v2, v3, v4, v5, v6)
    end
    function Vec6(X::Vec6)
        return new(X.v1, X.v2, X.v3, X.v4, X.v5, X.v6)
    end
    function Vec6(X::Array{<:Real,1})
        @assert length(X)==6
        return new(X...)
    end
end


macro Vec6(expr)
    @assert expr.head==:vect
    @assert length(expr.args)==6
    return esc(Expr(:call, :Vec6, expr.args...))
end

Base.convert(::Type{Vec6}, A::AbstractArray) = Vec6(A)


Base.length(::Vec6) = 6
Base.size(::Vec6) = return (6,)
Base.strides(::Vec6) = return (1,)
Base.unsafe_convert(::Type{Ptr{Float64}}, a::Vec6) = Base.unsafe_convert(Ptr{Float64}, pointer_from_objref(a))
Base.zeros(::Type{Vec6}) = return Vec6(0.0,0.0,0.0,0.0,0.0,0.0)
Base.copy(X::Vec6) = return Vec6(X.v1, X.v2, X.v3, X.v4, X.v5, X.v6)


function Base.copyto!(D::Vec6, S::Vec6)
    D.v1 = S.v1
    D.v2 = S.v2
    D.v3 = S.v3
    D.v4 = S.v4
    D.v5 = S.v5
    D.v6 = S.v6
    return D
end


function Base.getindex(X::Vec6, i::Int)
    if i<4
        i==1 && return X.v1
        i==2 && return X.v2
        return X.v3
    else
        i==4 && return X.v4
        i==5 && return X.v5
        return X.v6
    end

    throw(BoundsError(X,i))
end


function Base.setindex!(X::Vec6, val::Float64, i::Int)
    if i<4
        i==1 && (X.v1 = val)
        i==2 && (X.v2 = val)
        X.v3 = val
    else
        i==4 && (X.v4 = val)
        i==5 && (X.v5 = val)
        X.v6 = val
    end

    return val
end


Base.:*(X::Vec6, a::Real) = Vec6( a*X.v1, a*X.v2, a*X.v3, a*X.v4, a*X.v5, a*X.v6 )
Base.:*(a::Real, X::Vec6) = X*a
Base.:/(X::Vec6, a::Real) = X*(1.0/a)

Base.:*(M::Array{Float64,2}, X::Vec6) = Vec6( 
    M[1,1]*X.v1 + M[1,2]*X.v2 + M[1,3]*X.v3 + M[1,4]*X.v4 + M[1,5]*X.v5 + M[1,6]*X.v6,
    M[2,1]*X.v1 + M[2,2]*X.v2 + M[2,3]*X.v3 + M[2,4]*X.v4 + M[2,5]*X.v5 + M[2,6]*X.v6,
    M[3,1]*X.v1 + M[3,2]*X.v2 + M[3,3]*X.v3 + M[3,4]*X.v4 + M[3,5]*X.v5 + M[3,6]*X.v6,
    M[4,1]*X.v1 + M[4,2]*X.v2 + M[4,3]*X.v3 + M[4,4]*X.v4 + M[4,5]*X.v5 + M[4,6]*X.v6,
    M[5,1]*X.v1 + M[5,2]*X.v2 + M[5,3]*X.v3 + M[5,4]*X.v4 + M[5,5]*X.v5 + M[5,6]*X.v6,
    M[6,1]*X.v1 + M[6,2]*X.v2 + M[6,3]*X.v3 + M[6,4]*X.v4 + M[6,5]*X.v5 + M[6,6]*X.v6,
)

Base.:-(X::Vec6) = Vec6( -X.v1, -X.v2, -X.v3, -X.v4, -X.v5, -X.v6 )
Base.:+(X::Vec6, Y::Vec6) = Vec6( X.v1+Y.v1, X.v2+Y.v2, X.v3+Y.v3, X.v4+Y.v4, X.v5+Y.v5, X.v6+Y.v6 )
Base.:-(X::Vec6, Y::Vec6) = Vec6( X.v1-Y.v1, X.v2-Y.v2, X.v3-Y.v3, X.v4-Y.v4, X.v5-Y.v5, X.v6-Y.v6 )
Base.:+(X::Vec6, Y::Array{Float64,1}) = Vec6( X.v1+Y[1], X.v2+Y[2], X.v3+Y[3], X.v4+Y[4], X.v5+Y[5], X.v6+Y[6] )
Base.:-(X::Vec6, Y::Array{Float64,1}) = Vec6( X.v1-Y[1], X.v2-Y[2], X.v3-Y[3], X.v4-Y[4], X.v5-Y[5], X.v6-Y[6] )
Base.:+(X::Array{Float64,1}, Y::Vec6) = Vec6( X[1]+Y.v1, X[2]+Y.v2, X[3]+Y.v3, X[4]+Y.v4, X[5]+Y.v5, X[6]+Y.v6 )
Base.:-(X::Array{Float64,1}, Y::Vec6) = Vec6( X[1]-Y.v1, X[2]-Y.v2, X[3]-Y.v3, X[4]-Y.v4, X[5]-Y.v5, X[6]-Y.v6 )