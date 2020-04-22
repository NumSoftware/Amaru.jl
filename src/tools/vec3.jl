mutable struct Vec3<:AbstractArray{Float64,1}
    x::Float64
    y::Float64
    z::Float64
    function Vec3(x::Real=0.0, y::Real=0.0, z::Real=0.0)
        return new(x, y, z)
    end
    function Vec3(X::AbstractArray{<:Float64})
        return new(X[1], X[2], X[3])
    end
end

function Base.convert(::Type{Array{Float64,1}}, V::Vec3)
    return [ V.x, V.y, V.z]
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
