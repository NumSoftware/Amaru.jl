export Point, Line, Loop, Surface, GeoModel
export addpoint!, addline!, addarc!, addloop!, addsurface!, addvolume!

abstract type GeoEntity
end


mutable struct Point<:GeoEntity
    coord::Vec3
    adj::Array{Point,1}
    size::Float64
    id::Int
    tag::String
    function Point(x::Real, y::Real=0.0, z::Real=0.0; size=0.0, id=0, tag="")
        x = round(x, digits=8)
        y = round(y, digits=8)
        z = round(z, digits=8)
        return new(Vec3(x,y,z), Point[], size, id, tag)
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

Base.show(io::IO, obj::Point) = _show(io, obj, 2, "")

# Base.hash(n::Point) = hash( (round(n.coord.x+1.0, digits=8), round(n.coord.y+2.0, digits=8), round(n.coord.z+3.0, digits=8)) ) # see also hash(::Node)
Base.hash(n::Point) = UInt(n.id)
Base.:(==)(x::Point, y::Point) = (hash(x) == hash(y))
# Base.isequal(x::Point, y::Point) = (hash(x) == hash(y))


abstract type Curve<:GeoEntity
end

mutable struct Line<:Curve
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

mutable struct Arc<:Curve
    points::Array{Point,1} # second point is center
    surfaces::Array
    n::Int 
    id::Int
    tag::String

    function Arc(p1::Point, p2::Point, p3::Point; n=0, id=0)
        return new([p1, p2, p3], [], n, id)
    end
end



# Curves with same endpoints are considered equal
# thus, it is not possible to add two curves with the same endpoints

Base.hash(c::Curve) = hash( extrema((c.points[1].id, c.points[end].id)) )
Base.:(==)(x::Curve, y::Curve) = (hash(x) == hash(y))
Base.show(io::IO, c::Curve) = _show(io, c, 2, "")

# Base.hash(c::Curve) = hash( sort(collect(n.id for n in c.points[[1,end]])) )
# Base.hash(l::Curve) = hash( sort(collect(n.id for n in l.points)) )
# Base.:(==)(x::Line, y::Line) = (hash(x) == hash(y))
# Base.show(io::IO, line::Line) = _show(io, line, 2, "")


mutable struct Loop<:GeoEntity # related to CurveLoop
    curves::Array{Curve,1}
    id::Int

    function Loop(curves::Curve...; id=0)
        return new([curves...], id)
    end

    function Loop(curves::Array; id=0)
        return new(curves, id)
    end

end

Base.show(io::IO, loop::Loop) = _show(io, loop, 2, "")

mutable struct Surface<:GeoEntity
    loops::Array{Loop,1}
    volumes::Array
    id::Int
    tag::String
    transfinite::Bool
    recombine::Bool
    # id::Int

    function Surface(loops::Loop...; id=0, tag="")
        return new([loops...], Volume[], id, tag, false, false)
    end
end

Base.show(io::IO, surf::Surface) = _show(io, surf, 2, "")


mutable struct Volume<:GeoEntity # related to SurfaceLoop
    surfaces::Array{Surface,1}
    id::Int
    tag::String

    function Volume(surfs::Array{Surface,1}; tag="")
        return new(surfs, -1, tag)
    end
end

Base.show(io::IO, volume::Volume) = _show(io, volume, 2, "")




function getcoords(surf::Surface)
    return Float64[ p.coord[j] for p in surf.points, j=1:3 ]
end

function getcoords(points::Array{Point,1})
    return Float64[ p.coord[j] for p in points, j=1:3 ]
end


mutable struct GeoModel
    entities::OrderedSet{GeoEntity}
    _id::Int
    function GeoModel()
        return new(OrderedSet{GeoEntity}(), 0)
    end
end

Base.show(io::IO, geo::GeoModel) = _show(io, geo, 3, "")
function Base.getindex(s::OrderedSet, x) 
    x = getkey(s.dict, x, nothing)
    x === nothing && error("Key not found in OrderedSet")
    return x
end

function Base.getindex(set::OrderedSet{<:GeoEntity}, id::Int) 
    for ent in set
        ent.id == id && return ent
    end
    return nothing
end

function Base.getproperty(set::OrderedSet{<:GeoEntity}, s::Symbol)
    s==:points && return [ent for ent in set if ent isa Point]
    s==:lines && return [ent for ent in set if ent isa Line]
    s==:surfaces && return [ent for ent in set if ent isa Surface]
    s==:volumes && return [ent for ent in set if ent isa Volume]
    return getfield(set, s)
end


function getentity(geo::GeoModel, e::GeoEntity)
    return getkey(geo.entities.dict, e, nothing)
end

function getentity(geo::GeoModel, p::Point)
    pp = getkey(geo.entities.dict, p, nothing)
    pp===nothing || return pp
    
    tol = 1e-8
    for pp in geo.entities
        pp isa Point || continue
        norm(p.coord-pp.coord)<tol && return pp
    end

    return nothing
end

function getpoint(geo::GeoModel, p::Point)
    pp = getkey(geo.entities.dict, p, nothing)
    pp===nothing || return pp
    
    tol = 1e-8
    for pp in geo.entities
        pp isa Point || continue
        norm(p.coord-pp.coord)<tol && return pp
    end

    return nothing
end



@enum(GeoFindStatus,
    OUTSIDE  = 0,
    INSIDE   = 1,
    ATPOINT  = 2,
    ATBORDER = 3,
)

function Base.copy(geo::GeoModel, p::Point; dx=0.0, dy=0.0, dz=0.0)
    pp = Point(p.coord .+ [dx, dy, dz]; size=p.size, tag=p.tag)
    pp = addpoint!(geo, pp)
    return pp
end

function point_in_segment(X, X1, X2)
    tol = 1e-8
    X1X = X-X1
    X1X2 = X2-X1
    norm(cross(X1X2, X1X)) > tol && return false # check colinearity
    dot1 = dot(X1X2, X1X)
    dot2 = dot(X1X2, X1X2)
    return tol < dot1 < dot2-tol
end

function point_in_line(p::Point, l::Line)

    p in l.points[[1,end]] && return ATPOINT

    # tol = 1e-8
    X = p.coord
    X1 = l.points[1].coord
    X2 = l.points[end].coord

    X1X = X-X1
    X1X2 = X2-X1
    norm(cross(X1X2, X1X)) > tol && return false # check colinearity
    dot1 = dot(X1X2, X1X)
    dot2 = dot(X1X2, X1X2)

    0 < dot1 < dot2 && return INSIDE
    tol < dot1 < dot2-tol && return INSIDE
    # (dot1<tol || dot2>dot1+tol) && 

    return OUTSIDE
    # return tol < dot1 < dot2-tol
end



function addpoint!(geo::GeoModel, p::Point)
    pp = getpoint(geo, p)
    pp===nothing || return pp

    # add point
    geo._id +=1
    p.id = geo._id
    push!(geo.entities, p)

    # check if point in inside any existing line
    for l in geo.entities
        l isa Line || continue
        p1 = l.points[1]
        p2 = l.points[end]
        if point_in_segment(p.coord, p1.coord, p2.coord)
            # delete old and create two new lines
            l1 = addsingleline!(geo, p1, p, n=div(l.n,2), tag=l.tag)
            l2 = addsingleline!(geo, p, p2, n=div(l.n,2), tag=l.tag)

            for s in l.surfaces
                # fix loops that include l
                for lo in s.loops
                    idx = findfirst(==(l), lo.curves)
                    idx===nothing && continue
                    idxp1 = idx % length(lo.curves) + 1
                    # update loop lines
                    newlines = [ l1, l2 ]
                    p2 in lo.curves[idxp1].points || reverse!(newlines)
                    lo.curves = [ lo.curves[1:idx-1]; newlines; lo.curves[idx+1:end] ]
                end

                # update surface in l1 and l2
                push!(l1.surfaces, s)
                push!(l2.surfaces, s)
            end

            l.surfaces = []
            delete!(geo, l)
            break
        end
    end

    return p
end

function addpoint!(geo::GeoModel, x, y, z; size=0.0, tag="")
    return addpoint!(geo, Point(x,y,z; size=size, tag=tag))
end

function addpoint!(geo::GeoModel, X::Array; size=0.0, tag="")
    X = Vec3(X)
    return addpoint!(geo, Point(X[1],X[2],X[3]; size=size, tag=tag))
end


function coplanar(l1::Line, l2::Line)
    tol = 1e-8
    X1 = l1.points[1].coord
    V1 = l1.points[end].coord - X1

    X2 = l2.points[1].coord
    V2 = l2.points[end].coord - X2

    return abs(dot(cross(V1,V2), X1-X2)) < tol 
end

function intersection(l1::Line, l2::Line)
    # Assumes:
    # L1: X = X1 +V1*t
    # L2: Y = Y2 +V2*s
    # V1 = X2-X1
    # V2 = Y2-Y1

    tol = 1e-8
    X1 = l1.points[1].coord
    X2 = l1.points[end].coord
    V1 = X2-X1

    Y1 = l2.points[1].coord
    Y2 = l2.points[end].coord
    V2 = Y2-Y1

    #check if they are parallel or colinear
    norm(cross(V1, V2)) > tol || return nothing

    # check if lines are coplanar
    abs(dot(cross(V1,V2), X1-Y1)) < tol || return nothing

    # intersection point
    v1² = dot(V1,V1)
    v2² = dot(V2,V2)
    v1v2 = dot(V1,V2)
    s = ( dot(Y1-X1, V2)*v1² - dot(Y1-X1,V1)*v1v2) / (v1v2^2 - v1²*v2²)
    t = ( dot(Y1-X1, V1) + v1v2*s) / v1²

    # check if the intersection point is inside segments
    (-tol < s < 1+tol && -tol < t < 1+tol) || return nothing

    X = X1 + t*V1

    return Point(X...)
end

    
function findloops(p1, p2)
    loops = Loop[]
    (length(p1.adj)==1 || length(p2.adj)==1) && return loops
    
    visited = [ p1 ]
    prev = Dict()
    queue = [ p1 ]
    loops = []
    closing_segs = [  ]

    while length(queue)>0
        pp = popfirst!(queue)
        for p in pp.adj
            (haskey(prev, pp) && p == prev[pp]) && continue
            (pp,p) in closing_segs && continue
            if p in visited  # there is loop
                px = p
                loop = [ px ]
                while px!=p1
                    px = prev[px]
                    pushfirst!(loop, px)
                end

                px = pp
                push!(loop, pp)
                while px!=p1
                    px = prev[px]
                    push!(loop, px)
                end
                pop!(loop) # p1 was added twice

                p2 in loop || continue
                
                push!(loops, loop)
                push!(closing_segs, (p,pp))
                continue
            end
            push!(queue, p)
            push!(visited, p)
            prev[p] = pp
        end 
    end

    # @show length(loops)
    
    if length(loops)>2  # TODO: FIXME!
        loops = loops[1:2]
    end

    return loops
end


# function getcurveswithendpoint(geo::GeoModel, p::Point)
#     curves = []
#     for c in geo.entities
#         c isa Curve || continue
#         if p==c.points[1] || p==c.points[2]
#             push!(curves, p)
#         end
#     end
#     return curves
# end

# function getcurvefromendpoints(geo::GeoModel, p1::Point, p2::Point)
#     for c in geo.entities
#         c isa Curve || continue
#         endpoints = ( c.points[1], c.points[2] )
#         p1 in endpoints || continue
#         p2 in endpoints || continue
#         return c
#     end
#     return nothing
# end

function getcurve(geo::GeoModel, p1::Point, p2::Point)
    # all curves hatch function only considers endpoints
    # this function can return any curve type with the given endpoints
    l = Line(p1, p2)
    return getkey(geo.entities.dict, l, nothing)
end


function findloops(geo::GeoModel, curve::Curve)
    loops = Loop[]
    p1 = curve.points[1]
    p2 = curve.points[end]

    
    # @show p1.id
    # @show p2.id
    (length(p1.adj)==1 || length(p2.adj)==1) && return loops
    
    visited = [ p1 ]
    prev = Dict()
    queue = [ p1 ]
    loops = []
    visitededges = []

    function getsequence(p)
        seq = []
        p==p1 && return seq
        while true
            pr = prev[p]
            c = getcurve(geo, p, pr)
            push!(seq, c)
            p = pr
            pr==p1 && break
        end
        return seq
    end

    # put p2 as first adjacent of p1
    idxp2 = findfirst(x->x==p2, p1.adj)
    p1.adj = circshift(p1.adj, 1-idxp2)

    while length(queue)>0
        pp = popfirst!(queue) # pivot point
        # @show pp.id

        for p in pp.adj
            (haskey(prev, pp) && p==prev[pp]) && continue # do not look in previous points

            # advance while possible
            k=1
            pl = pp # last point in this path (pl)
            while true 
                p in visited && break
                length(p.adj)!=2 && break
                pnext = p.adj[1]==pl ? p.adj[2] : p.adj[1]
                prev[p] = pl
                push!(visited, p)
                pl = p
                p = pnext
                k = k+1
                # k==10 && error()
            end

            (pl,p) in visitededges && continue # skip if a loop was already closed with this segment

            push!(visitededges, (p,pl))
            push!(visitededges, (pl,p))

            # @show pl.id, p.id

            if p in visited  # there is loop
                # @show prev[p].id

                seq1 = reverse(getsequence(p))
                # @show 10
                seq2 = getsequence(pl)
                # @show 20

                c = getcurve(geo, p, pl)
                seq = [ seq1; c; seq2 ]
                # @show seq2

                # @show [ c.id for c in seq ]


                loop = Loop(seq)
                curve in seq || continue
                allunique(seq) || continue


                # @show [ c.id for c in loop.curves ]
                # @show [ p.id for p in seq1 ]
                # @show [ p.id for p in seq2 ]
                # @show [ p.id for p in getpoints(loop) ]

                # @show "yes"

                # set previous if seq2 contains initial curve
                if curve in seq2
                    # prev[p] = pl #??? seems to give problems (repeated segments)
                end

                push!(loops, loop)
                # push!(closing_segs, (p,pl))
                # push!(closing_segs, (pl,p))
                continue
            end
            push!(queue, p)
            push!(visited, p)
            prev[p] = pl
        end 
    end


    # check if there are loops with all points contained in other loops
    nloops = length(loops)
    indices = trues(nloops)
    for i in 1:nloops
        pts = getpoints(loops[i])
        for j in 1:nloops
            i==j && continue
            if issubset(pts, getpoints(loops[j]))
                indices[j] = false
                break
            end
        end
    end
    loops = loops[indices]

    # get only coplanar loops
    loops = [ lo for lo in loops if coplanar(lo) ]

    # check if resulting loops overlap each other
    nloops = length(loops)
    indices = trues(nloops)
    for i in 1:nloops
        for j in 1:nloops
            i==j && continue
            if overlaps(loops[i], loops[j])
                indices[j] = false
                break
            end
        end
    end
    loops = loops[indices]

    if length(loops)>2  # TODO: FIXME!
        loops = loops[1:2]
    end

    return loops
end


function getline(geo::GeoModel, p1::Point, p2::Point)
    l = Line(p1, p2)
    return getkey(geo.entities.dict, l, nothing)
end

# This fuction just adds a line without searching for loops
function addsingleline!(geo::GeoModel, p1::Point, p2::Point; n=0, tag="")
    ll = getline(geo, p1, p2)

    if ll!==nothing
        ll.n = n
        return ll
    end

    l = Line(p1, p2, tag=tag)
    geo._id +=1
    l.id = geo._id
    push!(geo.entities, l)
    push!(p1.adj, p2)
    push!(p2.adj, p1)

    return l
end

function addline!(geo::GeoModel, X1, X2; n=0, tag="")
    p1 = addpoint!(geo, X1)
    p2 = addpoint!(geo, X2)
    return addline!(geo, p1, p2; n=n, tag=tag)
end

function addline!(geo::GeoModel, p1::Point, p2::Point; n=0, tag="")
    l = Line(p1, p2)
    ll = getentity(geo, l)
    ll===nothing || return ll

    # look for intersections
    points = Point[ ]

    # check if line intersects other lines
    for li in geo.entities
        li isa Line || continue
        p = intersection(li, l)
        p === nothing && continue
        p = addpoint!(geo, p)
        push!(points, p)
    end

    # sort points    
    push!(points, p1)
    push!(points, p2)
    unique!(points)
    sort!(points, by=x->norm(x.coord-p1.coord))

    # add points
    for p in points
        addpoint!(geo, p) # may update some loops
    end

    # add lines
    npts = length(points)
    for i in 1:npts-1
        p1 = points[i]
        p2 = points[i+1]

        # l = Line(p1, p2, n=0, tag=tag)
        l = addsingleline!(geo, p1, p2, n=0, tag=tag)
        addloops!(geo, l)
    end
end

function addarc!(geo::GeoModel, p1::Point, p2::Point, p3::Point; n=0, tag="")
    p1 = getpoint(geo, p1)
    p2 = getpoint(geo, p2)
    p3 = getpoint(geo, p3)

    a = Arc(p1, p2, p3, n=n)
    aa = getentity(geo, a)
    aa===nothing || return aa

    # TODO: look for intersections

    # add arc
    geo._id +=1
    a.id = geo._id
    push!(geo.entities, a)
    push!(p1.adj, p3)
    push!(p3.adj, p1)
    addloops!(geo, a)

    return a
end

function addloops!(geo::GeoModel, c::Curve)
    loops = findloops(geo, c)
    for lo in loops
        addloop!(geo, lo) # may add a new surface
    end
end

function addline!2(geo::GeoModel, p1::Point, p2::Point; n=0, tag="")
    ll = getline(geo, p1, p2)

    ll===nothing || return ll

    l = Line(p1, p2)
    points = Point[ ]

    # check if line intersects other lines
    for li in geo.entities
        li isa Line || continue
        p = intersection(li, l)
        p === nothing && continue
        p = addpoint!(geo, p)
        push!(points, p)
    end

    # sort points    
    push!(points, p1)
    push!(points, p2)
    unique!(points)
    sort!(points, by=x->norm(x.coord-p1.coord))

    # add points
    for p in points
        addpoint!(geo, p) # may update some loops
    end

    # add lines
    npts = length(points)
    for i in 1:npts-1
        p1 = points[i]
        p2 = points[i+1]

        l = Line(p1, p2, n=0, tag=tag)
        ll = getkey(geo.entities.dict, l, nothing)
        ll===nothing || continue

        # add single lines
        geo._id +=1
        l.id = geo._id
        push!(geo.entities, l)
        push!(p1.adj, p2)
        push!(p2.adj, p1)

        # find coplanar loops
        loops = findloops(geo, l)
        loops = [ lo for lo in loops if coplanar(lo) ]

        # check if resulting loop overlap each other
        nloops = length(loops)
        indices = trues(nloops)
        for i in 1:nloops
            for j in 1:nloops
                i==j && continue
                if overlaps(loops[i], loops[j])
                    indices[j] = false
                    break
                end
            end
        end
        loops = loops[indices]

        # check if there two loops that overlaps the same surface
        if length(loops)==2
            ovs = nothing # overlapped surface
            for s in geo.entities
                s isa Surface || continue
                if overlaps(loops[1], s) && overlaps(loops[2], s)
                    ovs = s
                    break
                end
            end

            if ovs!==nothing
                # add loops manually
                for lo in loops
                    geo._id +=1
                    lo.id = geo._id
                    push!(geo.entities, lo)
                end

                s1 = addsurface!(geo, loops[1], tag=ovs.tag)
                s2 = addsurface!(geo, loops[2], tag=ovs.tag)
                for vol in ovs.volumes
                    push!(s1.volumes, vol)
                    push!(s2.volumes, vol)
                    push!(vol.surfaces, s1)
                    push!(vol.surfaces, s2)
                    idx = findfirst(s->s==ovs, vol.surfaces)
                    deleteat!(vol.surfaces, idx)
                end
                delete!(geo, ovs)

                continue
            end

        end

        for lo in loops
            addloop!(geo, lo) # may add a new surface
        end

    end

    return l
end


function Base.delete!(geo::GeoModel, l::Line)
    p1 = l.points[1]
    p2 = l.points[end]

    p1.adj = filter(x->x!=p2, p1.adj)
    p2.adj = filter(x->x!=p1, p2.adj)
    delete!(geo.entities, l)

    # delete loops associated with that line
    for s in l.surfaces
        for lo in s.loops
            l in lo.curves && delete!(geo, s)
        end
    end
end


function getarc(geo::GeoModel, p1::Point, p2::Point, p3::Point)
    a = Arc(p1, p2, p3)
    return getkey(geo.entities.dict, a, nothing)
end

# This fuction just adds an arc without searching for loops
function addsinglearc!(geo::GeoModel, p1::Point, p2::Point, p3::Point; n=0, tag="")
    aa = getarc(geo, p1, p2, p3)

    if aa!==nothing
        aa.n = n
        return aa
    end

    a = Arc(p1, p2, p3, tag=tag)
    geo._id +=1
    a.id = geo._id
    push!(geo.entities, a)
    push!(p1.adj, p3)
    push!(p3.adj, p1)

    return a
end


function getpoints(lo::Loop)
    points = Point[]
    l1 = lo.curves[1]

    # get the first point
    if lo.curves[1].points[end] in lo.curves[2].points[[1,end]]
        push!(points, l1.points[end])
    else
        push!(points, l1.points[1])
    end

    # get remaining points
    for i in 2:length(lo.curves)
        l = lo.curves[i]
        if l.points[1] == points[end]
            push!(points, l.points[end])
        else
            push!(points, l.points[1])
        end
    end

    return points
end

function coplanar(lo::Loop)
    tol = 1e-8
    points = getpoints(lo)
    
    # find a plane
    X1 = points[1].coord
    X2 = points[2].coord
    X1X2 = X2-X1

    # look for a non-colinear point
    local X, N
    for p in points[2:end]
        X = p.coord
        X1X = X-X1
        N = cross(X1X2, X1X)
        norm(N) > tol && break
    end
    
    # test the plane at each point
    for p in points
        X = p.coord

        if dot(X-X1, N) > tol
            return false
        end
    end

    return true
end

# loop points have to be coplanar
function getplane(lo::Loop)
    tol = 1e-8
    points = getpoints(lo)
    
    # find a plane
    X1 = points[1].coord
    X2 = points[2].coord
    X1X2 = X2-X1

    # look for a non-colinear point
    local N
    for p in points[2:end]
        X = p.coord
        X1X = X-X1
        N = cross(X1X2, X1X)
        norm(N) > tol && break
    end

    N = round.(N, digits=8)

    h = dot(X1,N)/dot(N,N)

    if h<0
        h = -h
        N = -N
    end

    return N, h
end

# checks if point is inside loop
function inside(p::Point, lo::Loop; withborder=true)
    N, h = getplane(lo)
    return inside(p, getpoints(lo), N; withborder)
end

function inside(p::Point, points::Array{Point,1}, normal::Array{Float64,1}; withborder=true)
    tol = 1e-8
    n = length(points)
    θ = 0.0
    normal = normalize(normal)
    for i in 1:n
        j = i*sign(n-i) + 1
        V1 = points[i].coord - p.coord
        V2 = points[j].coord - p.coord
        
        len1 = norm(V1)
        len2 = norm(V2)

        if len1*len2==0 # on vertex
            return withborder
        end

        val = clamp(dot(V1,V2)/(len1*len2), -1.0, 1.0) 

        N = normalize(cross(V1, V2))
        revsign = 1.0
        if norm(N-normal)>tol
            revsign = -1.0
        end

        # @show points[i].id
        # @show acos(val)
        θ += revsign*acos(val)
    end

    return abs(2*pi-abs(θ))<tol

end


# check if loop geometry overlaps another loop
function overlaps(lo1::Loop, lo2::Loop; withborder=true)
    tol = 1e-8
    N1, h1 = getplane(lo1)
    N2, h2 = getplane(lo2)

    # check if loop and surface are coplanar
    (h1!=h2 || norm(cross(N1,N2))>tol) && return false

    points = getpoints(lo1)

    # check if all loop points are inside the other loop
    for p in points
        inside(p, lo2, withborder=withborder) || return false
    end

    return true
end

# check if loop geometry overlaps a surface
function overlaps(lo::Loop, s::Surface; withborder=true)
    tol = 1e-8
    N1, h1 = getplane(lo)
    N2, h2 = getplane(s.loops[1])

    # check if loop and surface are coplanar
    (h1!=h2 || norm(cross(N1,N2))>tol) && return false

    points = getpoints(lo)
    # spoints = getpoints(s.loops[1])

    # check if all loop points are inside surface main loop
    for p in points
        inside(p, s.loops[1], withborder=withborder) || return false
        # inside(p, spoints, withborder=withborder) || return false

        # chek if p is outside hole loops
        if length(s.loops)>1
            for loo in s.loops[2:end]
                inside(p, loo, withborder=!withborder) && return false
            end
        end
    end

    return true
end


function addloop!(geo::GeoModel, lines::Line...)
    lo = Loop(lines...)
    return addloop!(geo, lo)
end

# function addentity(geo::GeoModel, ent)
#     found = getkey(geo.entities.dict, ent, nothing)
#     found===nothing || return ent

#     geo._id +=1
#     ent.id = geo._id
#     push!(geo.entities, ent)

#     return ent
# end

function addloop!(geo::GeoModel, lo::Loop)
    loo = getkey(geo.entities.dict, lo, nothing)
    loo===nothing || return loo

    geo._id +=1
    lo.id = geo._id
    push!(geo.entities, lo)

    # find if loop overlaps existing surface
    for s in geo.entities
        s isa Surface || continue

        # check if lo is a hole
        if overlaps(lo, s, withborder=false)
            # @show "hole"
            push!(s.loops, lo)
            for l in lo.curves
                push!(l.surfaces, s)
            end
            break
        end

        # check if lo is a subregion
        if overlaps(lo, s)
            # @show "ov"
            delete!(geo, s) # also deletes the loop
            break
        end
    end

    s = addsurface!(geo, lo)

    # check if there are loops that represet holes for this surface
    for loop in geo.entities
        loop isa Loop || continue

        if overlaps(loop, s, withborder=false)
            push!(s.loops, loop)
        end
    end

    return lo
end

function Base.delete!(geo::GeoModel, loop::Loop)
    delete!(geo.entities, loop)

    # delete surfaces associated with that loop
    for surf in geo.entities
        surf isa Surface || continue
        if length(surf.loops)==1
            loop==surf.loops[1] && delete!(geo, surf)
            continue
        else
            if loop==surf.loops[1] # remove external surface
                delete!(geo, surf)
            elseif loop in surf.loops[2:end] # remove holes
                surf.loops # TODO
            end
        end
    end

end


function addsurface!(geo::GeoModel, loops::Loop...; tag="")
    # @show "adding surface"
    s = Surface(loops..., tag=tag)
    ss = getkey(geo.entities.dict, s, nothing)
    ss === nothing || return ss

    geo._id +=1
    s.id = geo._id
    push!(geo.entities, s)

    # update edges
    for lo in s.loops
        for l in lo.curves
            push!(l.surfaces, s)
        end
    end

    return s
end

function Base.delete!(geo::GeoModel, surf::Surface)
    # @show "deleting surface"

    delete!(geo.entities, surf)

    for lo in surf.loops
        for l in lo.curves
            idx = findfirst(s->s==surf, l.surfaces)
            deleteat!(l.surfaces, idx)
        end

        # check if other surfaces uses this loop    
        found = false
        for s in geo.entities
            s isa Surface || continue
            if lo in s.loops
                found = true
                break
            end
        end
        found && continue
        delete!(geo.entities, lo)
    end

    # delete volumes associated
end


function addvolume!(geo::GeoModel, surfaces::Array{Surface,1}; tag="")
    v = Volume(surfaces, tag=tag)
    vv = getkey(geo.entities.dict, v, nothing)
    vv===nothing || return v

    geo._id +=1
    v.id = geo._id
    push!(geo.entities, v)

    return v
end

export extrude!


function extrude!(geo::GeoModel, line::Line; axis=[0,0,1], length=1.0)
    p1, p2 = line.points
    dx = length*axis[1]
    dy = length*axis[2]
    dz = length*axis[3]
    p4 = copy(geo, p1, dx=dx, dy=dy, dz=dz)
    p3 = copy(geo, p2, dx=dx, dy=dy, dz=dz)
    l1 = addsingleline!(geo, p1, p2, tag=line.tag)
    l2 = addsingleline!(geo, p2, p3, tag=line.tag)
    l3 = addsingleline!(geo, p3, p4, tag=line.tag)
    l4 = addsingleline!(geo, p4, p1, tag=line.tag)

    lo = Loop(l1, l2, l3, l4)

    geo._id +=1
    lo.id = geo._id
    push!(geo.entities, lo)

    s = addsurface!(geo, lo)
    return s
end

function extrude!(geo::GeoModel, surf::Surface; axis=[0,0,1], length=1.0)
    surfs = [ surf ]

    # extrude lateral lines
    for lo in surf.loops
        for line in lo.curves
            s = extrude!(geo, line, axis=axis, length=length)
            s.tag = surf.tag
            push!(surfs, s)
        end
    end

    # find lid loops
    loops = Loop[]
    for lo in surf.loops
        lines = Line[]
        for line in lo.curves
            points = Point[]
            for p in line.points
                pp = Point(p.coord .+ length.*axis)
                pp.tag = p.tag
                pp = addpoint!(geo, pp)
                push!(points, pp)
            end
            l = getline(geo, points...)
            push!(lines, l)
        end
        # lo = addloop!(geo, lines...)
        lo = Loop(lines...)
        push!(loops, lo)
    end

    for lo in loops
        geo._id +=1
        lo.id = geo._id
        push!(geo.entities, lo)
    end

    s = addsurface!(geo, loops..., tag=surf.tag)
    push!(surfs, s)

    v = addvolume!(geo, surfs, tag=surf.tag)

    for s in surfs
        push!(s.volumes, v)
    end

    return v
end

function extrude!(m::GeoModel, surfs::Array{Surface,1}; axis=[0,0,1], length=1.0)
    for s in surfs
        extrude!(m, s; axis=axis, length=length)
    end
end

export picksurface

function picksurface(geo::GeoModel, p::Point)
    return picksurface(geo, p.coord)
end

function picksurface(geo::GeoModel, coord)
    p = Point(coord)
    for s in geo.entities
        s isa Surface || continue
        isin = inside(p, s.loops[1])
        if isin && length(s.loops)>=2
            for lo in s.loops[2:end]
                if inside(p, lo)
                    isin = false
                    break
                end
            end
        end
        isin && return s
    end
    return nothing
end

function tag!(s::Surface, tag::String)
    s.tag = tag
end