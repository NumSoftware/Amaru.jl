mutable struct Vec3<:AbstractArray{Float64,1}
    x::Float64
    y::Float64
    z::Float64
    function Vec3(x::Float64, y::Float64, z::Float64)
        return new(x, y, z)
    end
    function Vec3(x::Real=0.0, y::Real=0.0, z::Real=0.0)
        return new(x, y, z)
    end
    function Vec3(X::Vec3)
        return new(X.x, X.y, X.z)
    end
    function Vec3(X::AbstractArray{<:Float64})
        return new(X[1], X[2], X[3])
    end
    function Vec3(X::Tuple{Float64,Float64,Float64})
        return new(X[1], X[2], X[3])
    end
end


const null_Vec3 = Vec3(Inf, Inf, Inf)
@inline null(::Type{Vec3}) = null_Vec3


function Base.convert(::Type{Array{Float64,1}}, V::Vec3)
    return [ V.x, V.y, V.z ]
end

function Base.convert(::Type{Vec3}, A::Array{Float64,1})
    n = length(A)
    @assert n in (1,2,3)
    if n==3
        return Vec3(A[1], A[2], A[3])
    elseif n==2
        return Vec3(A[1], A[2])
    else
        return Vec3(A[1])
    end
end

function Base.length(X::Vec3)
    return 3
end

function Base.size(X::Vec3)
    return (3,)
end

function Base.getindex(X::Vec3, i::Int)
    if i==1
        return X.x
    elseif i==2
        return X.y
    elseif i==3
        return X.z
    else
        throw(BoundsError(X,i))
    end
end

function Base.setindex!(X::Vec3, val, i::Int)
    if i==1
        return X.x = val
    elseif i==2
        return X.y = val
    elseif i==3
        return X.z = val
    else
        throw(BoundsError(X,i))
    end
end

function round!(X::Vec3; digits=0)
    X.x = round(X.x, digits=digits)
    X.y = round(X.y, digits=digits)
    X.z = round(X.z, digits=digits)
end

#macro commutative(expr)
    #expr2 = deepcopy(expr)
    #fargs = expr2.args[1].args
    #fargs[2], fargs[3] = fargs[3], fargs[2] 
    #return quote
        #$(esc(expr))
        #$(esc(expr2))
    #end
#end

Base.:*(X::Vec3, a::Real) = Vec3( X.x*a, X.y*a, X.z*a )
Base.:*(a::Real, X::Vec3) = Vec3( X.x*a, X.y*a, X.z*a )

Base.:+(X1::Vec3, X2::Vec3) = Vec3( X1.x+X2.x, X1.y+X2.y, X1.z+X2.z )
Base.:-(X1::Vec3, X2::Vec3) = Vec3( X1.x-X2.x, X1.y-X2.y, X1.z-X2.z )
Base.:+(X1::Vec3, X2::Array{Float64,1}) = Vec3( X1.x + X2[1], X1.y + X2[2], X1.z + X2[3] )
Base.:-(X1::Vec3, X2::Array{Float64,1}) = Vec3( X1.x - X2[1], X1.y - X2[2], X1.z - X2[3] )
Base.:+(X1::Array{Float64,1}, X2::Vec3) = Vec3( X1[1] + X2.x, X1[2] + X2.y, X1[3] + X2.z )
Base.:-(X1::Array{Float64,1}, X2::Vec3) = Vec3( X1[1] - X2.x, X1[2] - X2.y, X1[3] - X2.z )

function LinearAlgebra.normalize(X::Vec3)
    len = norm(X)
    return Vec3(X.x/len, X.y/len, X.z/len)
end
