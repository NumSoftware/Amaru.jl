export Point, Line, Loop, PlaneFace, Face, GeoModel
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

    function Point(coord::AbstractArray{<:Real}; size=0.0, id=0, tag="")
        x = round(coord[1], digits=8)
        y = round(coord[2], digits=8)
        z = length(coord)==2 ? 0.0 : round(coord[3], digits=8)
        return new(id, Vec3(x,y,z), GeoEntity[], size, tag)
    end
end


function Base.copy(p::Point)
    Point(p.coord, size=p.size, id=p.id)
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
    surfaces::AbstractArray
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
    surfaces::AbstractArray
    extrapoints::Array{Point,1} # points along the arc
    n::Int
    id::Int
    tag::String

    function Arc(p1::Point, p2::Point, p3::Point; n=0, id=0, tag="")
        this = new([p1, p2, p3], [], Point[], n, id, tag)

        # fill extra points along the arc; the second extra point is the middle
        v1 = p1.coord - p2.coord
        v2 = p3.coord - p2.coord
        θ = atan(norm(cross(v1,v2)), dot(v1,v2))
        α = θ/3
        â = normalize(v1)
        ĉ = normalize(cross(v1,v2))
        b̂ = normalize(cross(ĉ, â))
        r = norm(v1)
        for αi in [α, 0.5*θ, 2*α]
            p = p2.coord + r*cos(αi)*â + r*sin(αi)*b̂
            push!(this.extrapoints, Point(p))
        end

        return this
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

    function Loop(lines::AbstractArray; id=0)
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


abstract type AbstractFace<:GeoEntity
end


mutable struct PlaneFace<:AbstractFace
    id::Int
    loops::Vector{PlaneLoop}
    plane::Plane
    volumes::AbstractArray
    tag::String
    transfinite::Bool
    recombine::Bool

    function PlaneFace(loops::PlaneLoop...; id=0, tag="")
        return new(id, [loops...], loops[1].plane, Volume[], tag, false, false)
    end
end


mutable struct Face<:AbstractFace
    loops::Array{Loop,1}
    volumes::AbstractArray
    id::Int
    tag::String
    transfinite::Bool
    recombine::Bool

    function Face(loops::Loop...; id=0, tag="")
        return new([loops...], Volume[], id, tag, false, false)
    end
end


function Base.:(==)(s1::AbstractFace, s2::AbstractFace)
    return s1.loops[1]==s2.loops[1]
end

function getlines(s::AbstractFace)
    lines = AbstractLine[]
    for loop in s.loops
        for line in loop.lines
            push!(lines, line)
        end
    end
    return lines
end

function getpoints(s::AbstractFace)
    points = Set{Point}()
    for line in getlines(s)
        for point in line.points
            push!(points, point)
        end
    end
    return collect(points)
end


mutable struct Volume<:GeoEntity # related to SurfaceLoop
    surfaces::Array{AbstractFace,1}
    id::Int
    tag::String

    function Volume(surfs::Array{<:AbstractFace,1}; tag="")
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

function getpoints(v::Volume)
    points = Set{Point}()
    for surf in v.surfaces
        for point in getpoints(surf)
            push!(points, point)
        end
    end
    return collect(points)
end


mutable struct SurfaceLoop<:GeoEntity
    id::Int
    surfaces::Array{AbstractFace,1}

    function SurfaceLoop(surfaces::AbstractFace...; id=0)
        return new(id, [surfaces...])
    end
        
    function SurfaceLoop(surfaces::Vector{<:AbstractFace}; id=0)
        return SurfaceLoop(surfaces..., id=id)
    end
end



# Show functions
Base.show(io::IO, obj::Point) = _show(io, obj, 2)
Base.show(io::IO, c::AbstractLine) = _show(io, c, 2)
Base.show(io::IO, loop::Loop) = _show(io, loop, 2)
Base.show(io::IO, surf::AbstractFace) = _show(io, surf, 2)
Base.show(io::IO, volume::Volume) = _show(io, volume, 2)

