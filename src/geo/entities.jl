export Point, Line, CurvedLoop, FlatFace, CurvedFace, GeoModel
export addpoint!, addline!, addarc!, addloop!, addplanesurface!, addvolume!

abstract type GeoEntity
end

abstract type Edge<:GeoEntity
end


mutable struct Point<:GeoEntity
    id::Int
    coord::Vec3
    edges::Vector{<:Edge}
    size::Float64
    tag::String
    
    function Point(x::Real, y::Real=0.0, z::Real=0.0; size=0.0, id=0, tag="")
        x = round(x, digits=8)
        y = round(y, digits=8)
        z = round(z, digits=8)
        return new(id, Vec3(x,y,z), Edge[], size, tag)
    end

    function Point(coord::AbstractArray{<:Real}; size=0.0, id=0, tag="")
        x = round(coord[1], digits=8)
        y = round(coord[2], digits=8)
        z = length(coord)==2 ? 0.0 : round(coord[3], digits=8)
        return new(id, Vec3(x,y,z), Edge[], size, tag)
    end
end


function Base.copy(p::Point)
    Point(p.coord, size=p.size, id=0)
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



# Index operator for an collection of nodes
function Base.getindex(
    nodes::Array{Point,1}, 
    filters::Union{Expr,Symbolic,String}...;
    invert = false
    )
    
    idxs = collect(1:length(nodes))

    for filter in filters
        if typeof(filter) in (Expr, Symbolic)
            fnodes = nodes[idxs]

            T = Bool[]
            for node in fnodes
                x, y, z = node.coord.x, node.coord.y, node.coord.z
                push!(T, evaluate(filter, x=x, y=y, z=z))
            end
    
            idxs = idxs[T]
        elseif filter isa String
            idxs = [ i for i in idxs if nodes[i].tag==filter ]
        end
    end

    if invert
        idxs = setdiff(1:length(nodes), idxs)
    end

    return nodes[idxs]
end



struct Dart<:GeoEntity
    edge::Edge
    forward::Bool
    
    function Dart(edge::Edge; forward=true)
        return new(edge, forward)
    end
end



function Base.getproperty(dart::Dart, name::Symbol)
   if name==:points
        if dart.forward 
            return dart.edge.points
        else
            return reverse(dart.edge.points)
        end
    end

    getfield(dart, name)
end


function Base.:(==)(l1::Dart, l2::Dart) 
    return l1.forward == l2.forward && l1.edge == l2.edge
end


mutable struct Line<:Edge
    points::Vector{Point}
    faces::AbstractArray
    n::Int 
    id::Int
    tag::String
    
    function Line(points::Array{Point,1}; n=0, id=0, tag="")
        this = new(points, [], n, id, tag)
        return this
    end
end


function Line(p1::Point, p2::Point; n=0, id=0, tag="")
    return Line([p1, p2], n=n, id=id, tag=tag)
end


mutable struct Arc<:Edge
    id::Int
    points::Array{Point,1} # second point is center
    faces::AbstractArray
    extrapoints::Array{Point,1} # points along the arc
    n::Int
    tag::String

    function Arc(p1::Point, p2::Point, p3::Point; n=0, id=0, tag="")
        points = [p1, p2, p3]
        extrapoints = Point[]

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
            push!(extrapoints, Point(p))
        end

        this = new(id, points, Face[], extrapoints, n, tag)
        return this
    end
end


# Curves with same endpoints are considered equal
# thus, it is not possible to add two curves with the same endpoints
function Base.:(==)(l1::Edge, l2::Edge) 
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
    tol = 1e-6
    
    # find a plane
    X1 = points[1].coord
    X2 = points[2].coord
    X1X2 = X2-X1

    # look for a non-collinear point
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

    norm(N)<tol && error("Plane: Points are collinear")

    # test the plane at each point
    for p in points[k:end]
        X = p.coord
        if dot(X-X1, N) > tol
            @show dot(X-X1, N)
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


mutable struct Loop<:GeoEntity
    id::Int
    darts::Vector{Dart}
    flat::Bool

    function Loop(darts::Vector{Dart}; id=0, flat=false)
        return new(id, darts, flat)
    end
end


function Loop(edges::Vector{<:Edge}; id=0, flat=false)
    forward = edges[1].points[end] in edges[2].points[[1,end]]
    lastpoint = forward ? edges[1].points[end] : edges[1].points[1]
    dart1 = Dart(edges[1], forward=forward)

    darts = Dart[ dart1 ]
    for edge in edges[2:end]
        forward = lastpoint.id == edge.points[1].id
        lastpoint = forward ? edge.points[end] : edge.points[1]
        dart = Dart(edge, forward=forward)
        push!(darts, dart)
    end

    return Loop(darts, id=id, flat=flat)
end


function Loop(edges::Edge...; id=0, flat=false)
    return Loop([edges...]; id=id, flat=flat)
end


function Base.getproperty(loop::Loop, name::Symbol)
    if name==:edges
        return [ dart.edge for dart in loop.darts ]
    elseif name==:points
        points = Point[]
        for dart in loop.darts
            push!(points, dart.points[1])
        end
        return points
    end
  
    getfield(loop, name)
end





function Base.:(==)(lo1::Loop, lo2::Loop)
    length(lo1.darts)!=length(lo2.darts) && return false

    ids1 = [ l.edge.id for l in lo1.darts ]
    ids2 = [ l.edge.id for l in lo2.darts ]

    common = intersect(ids1, ids2)

    return length(common)==length(ids1)
end


mutable struct Face<:GeoEntity
    id::Int
    loops::Vector{Loop}
    volumes::AbstractArray
    flat::Bool
    tag::String
    transfinite::Bool
    recombine::Bool

    function Face(loops::Loop...; id=0, tag="")
        flat = loops[1].flat
        return new(id, [loops...], Volume[], flat, tag, false, false)
    end
end


function Base.:(==)(s1::Face, s2::Face)
    return s1.loops[1]==s2.loops[1]
end


struct FaceSpin
    face::Face
    normal::Vec3
    function FaceSpin(face::Face, normal::Vec3)
        return new(face, normal)
    end
end


function flip(spin::FaceSpin)
    return FaceSpin(spin.face, -spin.normal)
end


function Base.getproperty(face_spin::FaceSpin, name::Symbol)
    tol = 1e-6
    if name==:loops
        normal = getnormal(face_spin.face.loops[1])
        loops = Loop[]
        
        for (i,loop) in enumerate(face_spin.face.loops)
            # not reverse if collinear
            rev = norm(face_spin.normal - normal) < tol ? false : true
            
            # switch direction if hole
            i>1 && (rev = !rev)

            if rev
                push!(loops, reverse_loop(loop))
            else
                push!(loops, loop)
            end
        end
        return loops

     end
 
     getfield(face_spin, name)
 end
 
 
 function Base.:(==)(spin1::FaceSpin, spin2::FaceSpin) 
     return spin1.face == spin2.face && spin1.normal === spin2.normal
 end


mutable struct Volume<:GeoEntity # related to FaceLoop
    id::Int
    spins::Vector{FaceSpin}
    tag::String

    function Volume(spins::Vector{FaceSpin}; tag="", id=0)
        return new(id, spins, tag)
    end
end


function Base.:(==)(v1::Volume, v2::Volume)
    length(v1.spins)!=length(v2.spins) && return false
    ids1 = [ s.face.id for s in v1.spins ]
    ids2 = [ s.face.id for s in v2.spins ]
    common = intersect(ids1, ids2)
    return length(common)==length(ids1)
end


mutable struct FaceLoop<:GeoEntity
    id::Int
    spins::Vector{FaceSpin}

    function FaceLoop(spins::Vector{FaceSpin}; id=0)
        return new(id, spins)
    end
end




# Show functions
Base.show(io::IO, obj::Point) = _show(io, obj, 2)
Base.show(io::IO, edge::Edge) = _show(io, edge, 2)
Base.show(io::IO, loop::Loop) = _show(io, loop, 2)
Base.show(io::IO, surf::Face) = _show(io, surf, 2)
Base.show(io::IO, volume::Volume) = _show(io, volume, 2)
Base.show(io::IO, dart::Dart) = _show(io, dart, 2)


function Base.summary(dart::Dart) 
    str = "Dart: "
    str *= dart.forward ? "+" : "-"
    str *= "$(dart.edge.id)   "
    str *= "$(dart.points[1].id) -> $(dart.points[end].id)"
    str *= "   ·$(dart.edge.points[1].id)·$(dart.edge.points[end].id)"
end

function Base.show(io::IO, darts::Vector{Dart}) 
    print(io, "Vector{Dart}:  ")
    for dart in darts
        sgn = dart.forward ? "+" : "-"
        print(io, sgn, dart.edge.id)
    end
end

function Base.show(io::IO, spins::Vector{FaceSpin}) 
    print(io, "Vector{FaceSpin}:  ")
    for spin in spins
        print(io, " ", spin.face.id)
    end
end

function Base.summary(spin::FaceSpin) 
    str = "FaceSpin $(spin.face.id):"
    points = getpoints(spin.loops[1])
    for point in points
        str *= " $(point.id)"
    end
    return str
end

function Base.summary(loop::Loop) 
    str = "Loop $(loop.id):  "
    for dart in loop.darts
        str *= dart.forward ? "+" : "-"
        str *= "$(dart.edge.id)"
    end

    str *= "   "

    points = getpoints(loop)
    for point in points
        str *= "·$(point.id)"
    end
    return str
end