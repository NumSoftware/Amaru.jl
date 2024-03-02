export Point, Line, Loop, PlaneSurface, Surface, GeoModel
export addpoint!, addline!, addarc!, addloop!, addplanesurface!, addvolume!

abstract type GeoEntity
end

mutable struct Point<:GeoEntity
    id::Int
    coord::Vec3
    lines::Vector{<:GeoEntity}
    size::Float64
    tag::String
    function Point(x::Real, y::Real=0.0, z::Real=0.0; size=0.0, id=0, tag="")
        x = round(x, digits=8)
        y = round(y, digits=8)
        z = round(z, digits=8)
        return new(id, Vec3(x,y,z), GeoEntity[], size, tag)
    end
end


function Point(X::AbstractArray{<:Real}; size=0, id=0, tag="")
    if length(X)!=3
        X = [X; [0.0, 0.0]][1:3]
    end
    return Point(X...; size=size, id=id, tag=tag)
end


function Point(p::Point)
    Point(p.coord, p.size, id=p.id)
end


function Base.hash(n::Point) # required for unique function
    x = n.coord.x + 1
    y = n.coord.y + 2
    z = n.coord.z + 3
    xy = x*y
    xz = x*z
    yz = y*z
    return hash( (x, y, z, xy, xz, yz) )
end

Base.isequal(p1::Point, p2::Point) = hash(p1)==hash(p2)

function Base.:(==)(p1::Point, p2::Point)
    tol = 1e-8
    return norm(p1.coord-p2.coord)<tol
end


# Get node coordinates for a collection of nodes as a matrix
function getcoords(points::Array{Point,1}, ndim=3)
    npoints = length(points)
    return [ points[i].coord[j] for i in 1:npoints, j=1:ndim]
end


abstract type AbstractLine<:GeoEntity
end

mutable struct Line<:AbstractLine
    points::Array{Point,1}
    surfaces::Array
    n::Int 
    id::Int
    tag::String
    
    function Line(p1::Point, p2::Point; n=0, id=0, tag="")
        return new([p1, p2], [], n, id, tag)
    end

    function Line(points::Array{Point,1}; n=0, id=0, tag="")
        return new(points, [], n, id, tag)
    end
end

mutable struct Arc<:AbstractLine
    points::Array{Point,1} # second point is center
    surfaces::Array
    n::Int
    id::Int
    tag::String

    function Arc(p1::Point, p2::Point, p3::Point; n=0, id=0, tag="")
        return new([p1, p2, p3], [], n, id, tag)
    end
end



# Curves with same endpoints are considered equal
# thus, it is not possible to add two curves with the same endpoints
function Base.:(==)(l1::AbstractLine, l2::AbstractLine) 
    typeof(l1) != typeof(l2) && return false
    points1 = l1.points
    points2 = l2.points

    if points1[1].id!=points2[1].id 
        points2 = reverse(points2)
    end

    for i in 1:length(points1)
        points1[i].id!=points2[i].id && return false
    end
    
    return true
end



struct Plane
    normal::Vec3
    distance::Float64
end

function Base.:(==)(P1::Plane, P2::Plane) 
    tol = 1e-8
    norm(cross(P1.normal,P2.normal))<tol && abs(P1.distance-P2.distance)<tol && return true
end

Base.:(!=)(P1::Plane, P2::Plane) = !(P1==P2)

# create a plane
function Plane(N::Vec3, point::Point)
    N = round.(N, digits=8)
    N = normalize(N)
    h = dot(point.coord, N)

    if h<0
        h = -h
        N = -N
    end

    return Plane(N, h)
end

function Plane(points::Vector{Point})
    tol = 1e-8
    
    # find a plane
    X1 = points[1].coord
    X2 = points[2].coord
    X1X2 = X2-X1

    # look for a non-colinear point
    local N, k
    for i in 2:length(points)
        p = points[i]
        X = p.coord
        N = cross(X1X2, X-X1)
        if norm(N) > tol 
            k = i+1 # index to start to check if other points are in plane
            break
        end
    end

    norm(N)<tol && error("Plane: Points are colinear")

    # test the plane at each point
    for p in points[k:end]
        X = p.coord
        if dot(X-X1, N) > tol
            error("Plane: Points are not coplanar")
        end
    end

    N = round.(N, digits=8)
    N = normalize(N)
    h = dot(X1,N)

    if h<0
        h = -h
        N = -N
    end

    return Plane(N, h)
end


struct  DLine
    line::AbstractLine
    forward::Bool
end


abstract type AbstractLoop<:GeoEntity
end


mutable struct PlaneLoop<:AbstractLoop
    id::Int
    lines::Array{AbstractLine,1}
    plane::Plane

    function PlaneLoop(lines::AbstractLine...; id=0)
        points = [ p for l in lines for p in l.points ]
        plane = Plane(points)
        return new(id, [lines...], plane)
    end

    function PlaneLoop(lines::Vector{<:AbstractLine}; id=0)
        return PlaneLoop(lines..., id=id)
    end
end


mutable struct Loop<:AbstractLoop
    lines::Array{AbstractLine,1}
    id::Int

    function Loop(lines::AbstractLine...; id=0)
        return new([lines...], id)
    end

    function Loop(lines::Array; id=0)
        return new(lines, id)
    end
end


function Base.:(==)(lo1::AbstractLoop, lo2::AbstractLoop)
    length(lo1.lines)!=length(lo2.lines) && return false

    ids1 = [ l.id for l in lo1.lines ]
    ids2 = [ l.id for l in lo2.lines ]

    common = intersect(ids1, ids2)

    return length(common)==length(ids1)
end


abstract type AbstractSurface<:GeoEntity
end


mutable struct PlaneSurface<:AbstractSurface
    id::Int
    loops::Vector{PlaneLoop}
    plane::Plane
    volumes::Array
    tag::String
    transfinite::Bool
    recombine::Bool

    function PlaneSurface(loops::PlaneLoop...; id=0, tag="")
        return new(id, [loops...], loops[1].plane, Volume[], tag, false, false)
    end
end


mutable struct Surface<:AbstractSurface
    loops::Array{Loop,1}
    volumes::Array
    id::Int
    tag::String
    transfinite::Bool
    recombine::Bool

    function Surface(loops::Loop...; id=0, tag="")
        return new([loops...], Volume[], id, tag, false, false)
    end
end


function Base.:(==)(s1::AbstractSurface, s2::AbstractSurface)
    return s1.loops[1]==s2.loops[1]
end


mutable struct Volume<:GeoEntity # related to SurfaceLoop
    surfaces::Array{AbstractSurface,1}
    id::Int
    tag::String

    function Volume(surfs::Array{<:AbstractSurface,1}; tag="")
        return new(surfs, -1, tag)
    end
end


function Base.:(==)(v1::Volume, v2::Volume)
    length(v1.surfaces)!=length(v2.surfaces) && return false
    ids1 = [ s.id for s in v1.surfaces ]
    ids2 = [ s.id for s in v2.surfaces ]
    common = intersect(ids1, ids2)
    return length(common)==length(ids1)
end


# Show functions
Base.show(io::IO, obj::Point) = _show(io, obj, 2)
Base.show(io::IO, c::AbstractLine) = _show(io, c, 2)
Base.show(io::IO, loop::Loop) = _show(io, loop, 2)
Base.show(io::IO, surf::AbstractSurface) = _show(io, surf, 2)
Base.show(io::IO, volume::Volume) = _show(io, volume, 2)

