# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# Node
# =====
"""
A geometry type that represents a coordinate point.
"""
mutable struct Node
    x    ::Float64
    y    ::Float64
    z    ::Float64
    tag  ::String
    id   ::Int64
    _hash::UInt64
    function Node(x::Real, y::Real, z::Real=0.0; tag::String="")
        # zero is added to avoid negative bit sign for zero signed values
        x += 0.0
        y += 0.0
        z += 0.0

        _hash = hash( (round(x, digits=8), round(y, digits=8), round(z, digits=8)) )
        #x = round(x, digits=NDIG) + 0.0
        #y = round(y, digits=NDIG) + 0.0
        #z = round(z, digits=NDIG) + 0.0
        return new(x, y, z, tag,-1, _hash)
    end
    function Node(C::AbstractArray{<:Real}; tag::String="")
        # zero is added to avoid negative bit sign for zero signed values
        n = length(C)
        n==1 && return Node(C[1], 0.0, 0.0, tag=tag)
        n==2 && return Node(C[1], C[2], 0.0, tag=tag)
        return Node(C[1], C[2], C[3]+0.0, tag=tag)

        #n==1 && return new(C[1]+0.0, 0.0, 0.0, tag, -1)
        #n==2 && return new(C[1]+0.0, C[2]+0.0, 0.0, tag, -1)
        #return new(C[1]+0.0, C[2]+0.0, C[3]+0.0, tag, -1)
    end
end


### Node methods

Base.copy(point::Node) = Node(point.x, point.y, point.z, tag=point.tag)
Base.copy(points::Array{Node,1}) = [ copy(p) for p in points ]

# The functions below can be used in conjuntion with sort
get_x(point::Node) = point.x
get_y(point::Node) = point.y
get_z(point::Node) = point.z

#Base.hash(p::Node) = floor(UInt, abs(1000000 + p.coord.x*1001 + p.coord.y*10000001 + p.coord.z*100000000001))
#Base.hash(p::Node) = hash( (p.coord.x==0.0?0.0:p.coord.x, p.coord.y==0.0?0.0:p.coord.y, p.coord.z==0.0?0.0:p.coord.z) ) # comparisons used to avoid signed zero
#Base.hash(p::Node) = hash( (p.coord.x==0.0?0.0:p.coord.x, p.coord.y==0.0?0.0:p.coord.y, p.coord.z==0.0?0.0:p.coord.z) ) # comparisons used to avoid signed zero

Base.hash(p::Node) = hash( (round(p.coord.x, digits=8), round(p.coord.y, digits=8), round(p.coord.z, digits=8)) )
#Base.hash(p::Node) = hash( (p.coord.x, p.coord.y, p.coord.z))

Base.hash(points::Array{Node,1}) = sum(hash(p) for p in points)

Base.:(==)(p1::Node, p2::Node) = hash(p1)==hash(p2)

get_coords(p::Node) = [p.coord.x, p.coord.y, p.coord.z]

# Sentinel instance for not found point
const NO_POINT = Node(1e+308, 1e+308, 1e+308)

function get_node(points::Dict{UInt64,Node}, C::AbstractArray{<:Real})
    hs = hash(Node(C))
    return get(points, hs, nothing)
end

function get_coords(points::Array{Node,1}, ndim::Int64=3)
    np = length(points)
    C  = Array{Float64}(undef, np, ndim)
    for i=1:np
        C[i,1] = points[i].x
        C[i,2] = points[i].y
        ndim>2 && (C[i,3] = points[i].z)
    end

    return C
end

function setcoords!(points::Array{Node,1}, coords::AbstractArray{Float64,2})
    nrows, ncols = size(coords)
    @assert nrows == length(points)

    for (i,p) in enumerate(points)
        p.coord.x = coords[i,1]
        p.coord.y = coords[i,2]
        ncols==3 && (p.coord.z = coords[i,3])
    end
end


# Index operator for an collection of points
function Base.getindex(points::Array{Node,1}, filter_ex::Expr)
    R = Node[]
    for point in points
        x, y, z = point.x, point.y, point.z
        eval_arith_expr(filter_ex, x=x, y=y, z=z) && push!(R, point)
    end

    return R
end

function Base.getindex(points::Array{Node,1}, tag::String)
    return Node[ p for p in points if p.tag==tag ]
end

