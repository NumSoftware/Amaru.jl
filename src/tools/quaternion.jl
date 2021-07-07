mutable struct Quaternion<:AbstractArray{Float64,1}
    w::Float64
    x::Float64
    y::Float64
    z::Float64
    function Quaternion(w::Real, x::Real, y::Real, z::Real)
        return new(w, x, y, z)
    end
    function Quaternion(Q::AbstractArray{<:Float64})
        return new(Q[1], Q[2], Q[3], Q[4])
    end
    function Quaternion(Q::Tuple{Float64,Float64,Float64,Float64})
        return new(Q[1], Q[2], Q[3], Q[4])
    end
    function Quaternion(axis::AbstractArray{<:Float64}, angle::Real)
        return new(cos(angle/2), axis[1]*sin(angle/2), axis[2]*sin(angle/2), axis[3]*sin(angle/2))
    end
end


function Base.convert(::Type{Array{Float64,1}}, V::Quaternion)
    return [ V.w, V.x, V.y, V.z]
end

function Base.convert(::Type{Quaternion}, A::Array{Float64,1})
    n = length(A)
    @assert n==4
    return Quaternion(A[1], A[2], A[3], A[4])
end

function Base.length(Q::Quaternion)
    return 4
end

function Base.size(Q::Quaternion)
    return (4,)
end

function Base.getindex(Q::Quaternion, i::Int)
    if i==1
        return Q.w
    elseif i==2
        return Q.x
    elseif i==3
        return Q.y
    elseif i==4
        return Q.z
    else
        throw(BoundsError(Q,i))
    end
end

function Base.setindex!(Q::Quaternion, val, i::Int)
    if i==1
        return Q.w = val
    elseif i==2
        return Q.x = val
    elseif i==3
        return Q.y = val
    elseif i==4
        return Q.z = val
    else
        throw(BoundsError(Q,i))
    end
end

function Base.conj(Q::Quaternion)
    return Quaternion(Q.w, -Q.x, -Q.y, -Q.z)
end

function Base.:*(Q1::Quaternion, Q2::Quaternion)
    a1, b1, c1, d1 = Q1.w, Q1.x, Q1.y, Q1.z
    a2, b2, c2, d2 = Q2.w, Q2.x, Q2.y, Q2.z

    return Quaternion(
                      a1*a2 - b1*b2 - c1*c2 - d1*d2,
                      a1*b2 + b1*a2 + c1*d2 - d1*c2,
                      a1*c2 - b1*d2 + c1*a2 + d1*b2,
                      a1*d2 + b1*c2 - c1*b2 + d1*a2
                     )
end

function Base.:*(Q::Quaternion, V::Vec3)
    a1, b1, c1, d1 = Q.w, Q.x, Q.y, Q.z
    a2, b2, c2, d2 = 0.0, V.x, V.y, V.z

    return Quaternion(
                      a1*a2 - b1*b2 - c1*c2 - d1*d2,
                      a1*b2 + b1*a2 + c1*d2 - d1*c2,
                      a1*c2 - b1*d2 + c1*a2 + d1*b2,
                      a1*d2 + b1*c2 - c1*b2 + d1*a2
                     )
end

Base.:+(V::Vec3, Q::Quaternion) = Vec3( V.x+Q.x, V.y+Q.y, V.z+Q.z )
Base.:+(Q::Quaternion, V::Vec3) = Vec3( V.x+Q.x, V.y+Q.y, V.z+Q.z )
Base.:-(V::Vec3, Q::Quaternion) = Vec3( V.x-Q.x, V.y-Q.y, V.z-Q.z )
Base.:-(Q::Quaternion, V::Vec3) = Vec3( Q.x-V.x, Q.y-V.y, Q.z-V.z )

function rotate(X::Array{<:Real,1}, Q::Quaternion)
    Y = Q*X*conj(Q)
    return Vec3(Y[2], Y[3], Y[4])
end