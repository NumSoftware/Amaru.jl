export Point, Line, Loop, Surface, GeoModel
export addpoint!, addline!, addloop!, addsurface!, addvolume!

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




abstract type AbstractCurve<:GeoEntity
end

mutable struct Line<:AbstractCurve
    points::Array{Point,1}
    surfaces::Array
    n::Int 
    id::Int
    
    function Line(p1::Point, p2::Point; n=0, id=0)
        return new([p1, p2], [], n, id)
    end

    function Line(points::Array{Point,1}; n=0, id=0)
        return new(points, [], n, id)
    end
end

mutable struct Arc<:AbstractCurve
    points::Array{Point,1} # second point is center
    id::Int

    function Arc(p1::Point, pc::Point, p2::Point; n=0, id=0)
        return new([p1, pc, p2], id)
    end
end

# Base.hash(l::AbstractCurve) = 
# begin
#     sum( hash(p) for p in l.points)
# end

Base.hash(l::AbstractCurve) = hash( sort(collect(n.id for n in l.points)) )

Base.:(==)(x::Line, y::Line) = (hash(x) == hash(y))

Base.show(io::IO, line::Line) = _show(io, line, 2, "")


mutable struct Loop<:GeoEntity # related to CurveLoop
    lines::Array{AbstractCurve,1}
    id::Int

    function Loop(lines::AbstractCurve...; id=0)
        return new([lines...], id)
    end

end

Base.show(io::IO, loop::Loop) = _show(io, loop, 2, "")

mutable struct Surface<:GeoEntity
    loops::Array{Loop,1}
    volumes::Array
    id::Int

    function Surface(loops::Loop...; id=0)
        return new([loops...], Volume[], id)
    end
end

Base.show(io::IO, surf::Surface) = _show(io, surf, 2, "")



mutable struct Volume<:GeoEntity # related to SurfaceLoop
    surfaces::Array{Surface,1}
    id::Int

    function Volume(surfs::Array{Surface,1})
        return new(surfs)
    end

 end


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




@enum(GeoFindStatus,
    OUTSIDE  = 0,
    INSIDE   = 1,
    ATPOINT  = 2,
    ATBORDER = 3,
)

function Base.copy(m::GeoModel, p::Point; dx=0.0, dy=0.0, dz=0.0)
    pp = Point(p.coord .+ [dx, dy, dz]; size=p.size, tag=p.tag)
    pp = addpoint!(m, pp)
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


function addpoint!(geo::GeoModel, p::Point)
    pp = getpoint(geo, p)
    # pp = getkey(geo.entities.dict, p, nothing)
    pp===nothing || return pp

    # check if point is close to existing points
    # for pp in geo.entities
        # pp isa Point || continue
        # tol = 1e-8
        # norm(p.coord-pp.coord)<tol && return pp
    # end

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
            l1 = addlinenoloops(geo, p1, p, n=div(l.n,2))
            l2 = addlinenoloops(geo, p, p2, n=div(l.n,2))

            for s in l.surfaces
                # fix loops that include l
                for lo in s.loops
                    idx = findfirst(==(l), lo.lines)
                    idx===nothing && continue
                    idxp1 = idx % length(lo.lines) + 1
                    # update loop lines
                    newlines = [ l1, l2 ]
                    p2 in lo.lines[idxp1].points || reverse!(newlines)
                    lo.lines = [ lo.lines[1:idx-1]; newlines; lo.lines[idx+1:end] ]
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

function addpoint!(geo::GeoModel, x, y, z; size=0.0)
    return addpoint!(geo, Point(x,y,z; size=size))
end

function addpoint!(geo::GeoModel, X::Array; size=0.0)
    X = Vec3(X)
    return addpoint!(geo, Point(X[1],X[2],X[3]; size=size))
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

    return loops
end

function getline(geo::GeoModel, p1::Point, p2::Point)
    l = Line(p1, p2)
    return getkey(geo.entities.dict, l, nothing)
end

function addlinenoloops(geo::GeoModel, p1::Point, p2::Point; n=0)
    ll = getline(geo, p1, p2)
    # l = Line(p1, p2, n=n)
    # ll = getkey(geo.entities.dict, l, nothing)

    if ll!==nothing
        ll.n = n
        return ll
    end

    l = Line(p1, p2)
    geo._id +=1
    l.id = geo._id
    push!(geo.entities, l)
    push!(p1.adj, p2)
    push!(p2.adj, p1)

    return l
end

function addline!(geo::GeoModel, X1, X2; n=0)
    p1 = addpoint!(geo, X1)
    p2 = addpoint!(geo, X2)
    return addline!(geo, p1, p2; n=n)
end

function addline!(geo::GeoModel, p1::Point, p2::Point; n=0)
    ll = getline(geo, p1, p2)
    # l = Line(p1, p2, n=n)
    # ll = getkey(geo.entities.dict, l, nothing)

    if ll!==nothing
        ll.n = n
        return ll
    end

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

        l = Line(p1, p2, n=0)
        ll = getkey(geo.entities.dict, l, nothing)
    
        ll===nothing || continue

        geo._id +=1
        l.id = geo._id
        push!(geo.entities, l)
        push!(p1.adj, p2)
        push!(p2.adj, p1)

        # find loops pts if any loops
        loops_pts = findloops(p1, p2)
        loops = Loop[]

        # find plane loops
        
        for pts in loops_pts
            # convert pts to loop
            lines = AbstractCurve[]
            n = length(pts)
            for i in 1:n
                j = i*sign(n-i) + 1
                l = Line(pts[i], pts[j])
                l = getkey(geo.entities.dict, l, nothing)
                push!(lines, l)
            end

            # add loop
            lo = Loop(lines...)
            coplanar(lo) && push!(loops, lo)
        end

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

                s1 = addsurface!(geo, loops[1])
                s2 = addsurface!(geo, loops[2])
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
            addloop!(geo, lo)
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
            l in lo.lines && delete!(geo, s)
        end
    end
end


function getpoints(lo::Loop)
    points = Point[]
    l1 = lo.lines[1]

    # get the first point
    if lo.lines[1].points[end] in lo.lines[2].points[[1,end]]
        push!(points, l1.points[end])
    else
        push!(points, l1.points[1])
    end

    # get remaining points
    for i in 2:length(lo.lines)
        l = lo.lines[i]
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
    return inside(p, getpoints(lo); withborder)
end

function inside(p::Point, points::Array{Point,1}; withborder=true)
    tol = 1e-8
    n = length(points)
    θ = 0.0
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
        θ += acos(val)
    end

    return abs(2*pi-θ)<tol

end


# check if loop geometry overlaps a surface
function overlaps(lo::Loop, s::Surface; withborder=true)
    tol = 1e-8
    N1, h1 = getplane(lo)
    N2, h2 = getplane(s.loops[1])

    # check if loop and surface are coplanar
    (h1!=h2 || norm(cross(N1,N2))>tol) && return false

    points = getpoints(lo)
    spoints = getpoints(s.loops[1])

    # check if all loop points are inside surface main loop
    for p in points
        inside(p, spoints, withborder=withborder) || return false

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
            for l in lo.lines
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


function addsurface!(geo::GeoModel, loops::Loop...)
    # @show "adding surface"
    s = Surface(loops...)
    ss = getkey(geo.entities.dict, s, nothing)
    ss === nothing || return ss

    geo._id +=1
    s.id = geo._id
    push!(geo.entities, s)

    # update edges
    for lo in s.loops
        for l in lo.lines
            push!(l.surfaces, s)
        end
    end

    return s
end

function Base.delete!(geo::GeoModel, surf::Surface)
    # @show "deleting surface"

    delete!(geo.entities, surf)

    for lo in surf.loops
        for l in lo.lines
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


function addvolume!(geo::GeoModel, surfaces::Array{Surface,1})
    v = Volume(surfaces)
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
    # p4 = Point(p1.coord .+ length.*axis)
    # p3 = Point(p2.coord .+ length.*axis)
    # p4 = addpoint!(m, p4)
    # p3 = addpoint!(m, p3)
    p4 = copy(geo, p1, dx=dx, dy=dy, dz=dz)
    p3 = copy(geo, p2, dx=dx, dy=dy, dz=dz)
    l1 = addlinenoloops(geo, p1, p2)
    l2 = addlinenoloops(geo, p2, p3)
    l3 = addlinenoloops(geo, p3, p4)
    l4 = addlinenoloops(geo, p4, p1)

    # lo = addloop!(m, l1, l2, l3, l4)
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
        for line in lo.lines
            s = extrude!(geo, line, axis=axis, length=length)
            push!(surfs, s)
        end
    end

    # find lid loops
    loops = Loop[]
    for lo in surf.loops
        lines = Line[]
        for line in lo.lines
            points = Point[]
            for p in line.points
                pp = Point(p.coord .+ length.*axis)
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

    s = addsurface!(geo, loops...)
    push!(surfs, s)

    v = addvolume!(geo, surfs)

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
