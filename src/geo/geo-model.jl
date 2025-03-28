export Point, Line, Loop, PlaneFace, Face, GeoModel
export addpoint!, addline!, addpolygon!, addarc!, addloop!, addplanesurface!, addvolume!, addblock!

mutable struct GeoModel
    points::Vector{Point}
    lines::Vector{AbstractLine}
    loops::Vector{AbstractLoop}
    surfaces::Vector{AbstractFace}
    volumes::Vector{Volume}
    subpaths::Vector{SubPath}
    blocks::Vector{AbstractBlock}
    size::Float64
    ndim::Int
    _id::Int

    function GeoModel(; size=0.0)
        return new( [], [], [], [], [], [], [], size, 0, 0 )
    end
end

Base.show(io::IO, geo::GeoModel) = _show(io, geo, 3, "")


function Base.getindex(geo::GeoModel, idx::Int)
    for p in geo.points
        p.id==idx && return p
    end
    for l in geo.lines
        l.id==idx && return l
    end
    for lo in geo.loops
        lo.id==idx && return lo
    end
    for s in geo.surfaces
        s.id==idx && return s
    end
    for v in geo.volumes
        v.id==idx && return v
    end
    for sp in geo.subpaths
        sp.id==idx && return sp
    end
    for b in geo.blocks
        b.id==idx && return b
    end
    return nothing
end



function addblock!(geo::GeoModel, block::AbstractBlock)
    push!(geo.blocks, block)
end

function getpoint(geo::GeoModel, p::Point)
    for pp in geo.points
        p==pp && return pp
    end

    return nothing
end

function getline(geo::GeoModel, line::AbstractLine)
    idx = findfirst(==(line), geo.lines)
    idx === nothing && (idx=0) 
    return get(geo.lines, idx, nothing)
end

function getline(geo::GeoModel, p1::Point, p2::Point)
    l = Line(p1, p2)
    return getline(geo, l)
end

function getloop(geo::GeoModel, loop::AbstractLoop)
    idx = findfirst(==(loop), geo.loops)
    idx === nothing && return idx
    return geo.loops[idx]
end

function getsurface(geo::GeoModel, surf::AbstractFace)
    idx = findfirst(==(surf), geo.surfaces)
    idx === nothing && return idx
    return geo.surfaces[idx]
end

function getvolume(geo::GeoModel, vol::Volume)
    idx = findfirst(==(vol), geo.volumes)
    idx === nothing && return idx
    return geo.volumes[idx]
end



function join_volumes!(geo::GeoModel, v1::Volume, v2::Volume)

    # select find inner and outer surfaces
    innerfaces = intersect(v1.surfaces, v2.surfaces)
    outerfaces = setdiff(union(v1.surfaces, v2.surfaces), innerfaces)

    # remove innerfaces
    for s in innerfaces
        delete!(geo, s)
    end

    # add new volume
    v = addvolume!(geo, collect(outerfaces))

    return v
end


function join_planar_surfaces!(geo::GeoModel, s1::PlaneFace, s2::PlaneFace)
    # check if surfaces are coplanar
    s1.plane != s2.plane && return

    holes = [ s1.loops[2:end]; s2.loops[2:end] ]

    # get inner and outer lines
    innerlines = intersect(s1.loops[1].lines, s2.loops[1].lines)
    outerlines = setdiff(union(s1.loops[1].lines, s2.loops[1].lines), innerlines)

    # remove inner lines
    for l in innerlines
        delete!(geo, l)
    end

    # find loop
    loops = findloops(outerlines[1], lines=outerlines)
    s = addplanesurface!(geo, loops[1])

    return s
end


function Base.delete!(geo::GeoModel, line::Line)
    filter!(!=(line), geo.lines)

    p1 = line.points[1]
    p2 = line.points[end]

    filter!(!=(line), p1.lines)
    filter!(!=(line), p2.lines)

    # check if line is between two coplanar faces
    if length(line.surfaces)==2 && line.surfaces[1].plane==line.surfaces[2].plane
        s1 = line.surfaces[1] # to be deleted
        s2 = line.surfaces[2] # to be extended

        holes = [ s1.loops[2:end]; s2.loops[2:end] ]

        # get inner and outer lines
        innerlines = intersect(s1.loops[1].lines, s2.loops[1].lines)
        outerlines = setdiff(union(s1.loops[1].lines, s2.loops[1].lines), innerlines)

        # find loop and update s2 outer loop
        loops = findloops(outerlines[1], lines=outerlines)
        loop = addsingleloop!(geo, loops[1])
        s2.loops[1] = loop

        # remove s1 from geo
        filter!(!=(s1), geo.surfaces)
        filter!(!=(s1.loops[1]), geo.loops)

        # remove references to s1 in s1 lines and update with s2
        for lo in s1.loops
            for l in lo.lines
                l in innerlines && continue
                filter!(!=(s1), l.surfaces)
                push!(l.surfaces, s2)
            end
        end

        # remove references to s1 and s2 in inner lines
        for l in innerlines
            filter!(!=(s1), l.surfaces)
            filter!(!=(s2), l.surfaces)
        end

        # remove s1 from volumes associated to s2
        for v in s2.volumes
            filter!(!=(s1), v.surfaces)
        end

        # fix holes
    
    else
        # delete surfaces associated with that line
        for s in line.surfaces
            for lo in s.loops
                line in lo.lines && delete!(geo, s)
            end
        end
    end
end


function Base.delete!(geo::GeoModel, surf::AbstractFace)


    filter!(!=(surf), geo.surfaces)

    if !isa(surf, PlaneFace) 
        # delete loop if not plane surface
        delete!(geo.loops, lo) # todo: check
    else
        # remove references
        for lo in surf.loops
            # for each curve remove links to this surface
            for l in lo.lines
                filter!(!=(surf), l.surfaces)
            end 
    
            # check if other surfaces uses this loop
            found = false
            for s in geo.surfaces
                s.id==surf.id && continue
                s isa PlaneFace || continue
                if lo in s.loops
                    found = true
                    break
                end
            end
    
            if !found
                filter!(!=(lo), geo.loops)
            end
            # !found && delete!(geo.loops, lo) # todo: check
        end
    end

    # remove volumes that uses this surface

    
    if length(surf.volumes)==2
        v1 = surf.volumes[1] # volume to be removed
        v2 = surf.volumes[2] # volume to be extended

        innerfaces = intersect(v1.surfaces, v2.surfaces)
        outerfaces = setdiff(union(v1.surfaces, v2.surfaces), innerfaces)

        # update faces in v2
        v2.surfaces = outerfaces

        # remove v1 from geo
        filter!(!=(v1), geo.volumes)

        # remove references to v1 in v1 surfaces and update with v2
        for s in v1.surfaces
            s in innerfaces && continue
            filter!(!=(v1), s.volumes)
            push!(s.volumes, v2)
        end

        # remove v1 and v2 from innerfaces
        for s in innerfaces
            filter!(!=(v1), s.volumes)
            filter!(!=(v2), s.volumes)
        end

        # to_join && addvolume!(geo, collect(outerfaces))
    else
        for v in surf.volumes
            delete!(geo, v)
        end
        # vols = Volume[ v for v in geo.volumes if surf in v.surfaces ]
        # geo.volumes = setdiff(geo.volumes, vols)
    end

end

function Base.delete!(geo::GeoModel, v::Volume)
    filter!(!=(v), geo.volumes)

    # filter v from v surfaces
    for s in v.surfaces
        filter!(!=(v), s.volumes)
    end

end


function insegment(X, X1, X2)
    tol = 1e-8
    X1X = X-X1
    X1X2 = X2-X1
    norm(cross(X1X2, X1X)) > tol && return false # check colinearity
    dot1 = dot(X1X2, X1X)
    dot2 = dot(X1X2, X1X2)
    return tol < dot1 < dot2-tol # check if inside segment within tol
end

function insegment(p::Point, p1::Point, p2::Point)
    return insegment(p.coord, p1.coord, p2.coord)
end


@enum(GeoFindStatus,
    OUTSIDE  = 0,
    INSIDE   = 1,
    ATPOINT  = 2,
    ATBORDER = 3,
)

function point_in_line(p::Point, l::Line)

    p in l.points[[1,end]] && return ATPOINT

    tol = 1e-8
    X = p.coord
    X1 = l.points[1].coord
    X2 = l.points[end].coord

    X1X = X-X1
    X1X2 = X2-X1
    norm(cross(X1X2, X1X)) > tol && return false # check colinearity
    dot1 = dot(X1X2, X1X)
    dot2 = dot(X1X2, X1X2)

    0 < dot1 < dot2 && return INSIDE

    return OUTSIDE
end


function coplanar(P::Plane, l::AbstractLine)
    tol = 1e-6 # reduced tolerance 

    for p in l.points # check if every point is in the plane
        abs(dot(P.normal, p.coord) - P.distance) > tol && return false
    end

    return true
end


function coplanar(P::Plane, p::Point)
    tol = 1e-8 
    abs(dot(P.normal, p.coord) - P.distance) < tol && return true
    return false
end

function addsinglepoint!(geo::GeoModel, pt::Point)
    p = getpoint(geo, pt)
    p===nothing || return p
    
    geo._id +=1
    pt.id = geo._id
    push!(geo.points, pt)

    return pt
end


function addpoint!(geo::GeoModel, pt::Point)
    p = getpoint(geo, pt)
    p===nothing || return p

    # add point
    geo._id +=1
    pt.id = geo._id
    push!(geo.points, pt)

    geo.ndim = pt.coord[3] != 0.0 ? max(geo.ndim, 2) : 3

    # check if point is inside any existing line
    for l in geo.lines
        p1 = l.points[1]
        p2 = l.points[end]
        if insegment(pt.coord, p1.coord, p2.coord)
            newline = addsingleline!(geo, pt, p2, n=div(l.n,2), tag=l.tag) # create a new line
            l.points = [ p1, pt ] # updagte original line
            pt.lines = GeoEntity[ l, newline ] # update pt links
            filter!(!=(l), p2.lines) # remove ref to l in p2

            # update existing loops that contained the original line
            for lo in geo.loops
                idx = findfirst(==(l), lo.lines)
                idx===nothing && continue
                idxp1 = idx % length(lo.lines) + 1
                bothlines = [ l, newline ]
                p2 in lo.lines[idxp1].points || reverse!(bothlines)
                lo.lines = [ lo.lines[1:idx-1]; bothlines; lo.lines[idx+1:end] ]
            end

            newline.surfaces = copy(l.surfaces)
            break
        end
    end

    return pt
end


function addpoint!(geo::GeoModel, x, y, z=0.0; size=0.0, tag="")
    return addpoint!(geo, Point(x,y,z; size=size, tag=tag))
end


function addpoint!(geo::GeoModel, X::AbstractArray; size=0.0, tag="")
    X = Vec3(X)
    return addpoint!(geo, Point(X[1],X[2],X[3]; size=size, tag=tag))
end

function Base.copy!(geo::GeoModel, p::Point; dx=0.0, dy=0.0, dz=0.0)
    pp = Point(p.coord .+ [dx, dy, dz]; size=p.size, tag=p.tag)
    pp = addpoint!(geo, pp)
    return pp
end

# This fuction just adds a line without searching for loops
function addsingleline!(geo::GeoModel, p1::Point, p2::Point; n=0, tag="")
    l = getline(geo, p1, p2)

    if l !== nothing
        l.n = n
        return l
    end

    line = Line(p1, p2, n=n, tag=tag)

    geo._id +=1
    line.id = geo._id

    push!(geo.lines, line)
    push!(p1.lines, line)
    push!(p2.lines, line)
    return line
end


# This fuction just adds an arc without searching for loops
function addsinglearc!(geo::GeoModel, p1::Point, p2::Point, p3::Point; n=0, tag="")
    p1 = getpoint(geo, p1)
    p2 = getpoint(geo, p2)
    p3 = getpoint(geo, p3)

    arc = Arc(p1, p2, p3, n=n, tag=tag)
    a = getline(geo, arc)
    a===nothing || return a

    # add arc
    geo._id +=1
    arc.id = geo._id

    push!(geo.lines, arc)
    push!(p1.lines, arc)
    push!(p3.lines, arc)

    return arc
end


function addsingleloop!(geo::GeoModel, loop::AbstractLoop)
    lo = getloop(geo, loop)
    lo===nothing || return lo

    geo._id +=1
    loop.id = geo._id
    push!(geo.loops, loop)

    return loop
end


function intersection(l1::Line, l2::Line)
    # Assumes:
    # L1: P = P1 +V1*t
    # L2: Q = Y2 +V2*s
    # V1 = P2-P1
    # V2 = Q2-Q1

    tol = 1e-8
    P1 = l1.points[1].coord
    P2 = l1.points[end].coord
    V1 = P2-P1

    Q1 = l2.points[1].coord
    Q2 = l2.points[end].coord
    V2 = Q2-Q1

    #check if they are parallel or colinear
    norm(cross(V1, V2)) > tol || return nothing

    # check if lines are coplanar
    abs(dot(cross(V1,V2), P1-Q1)) < tol || return nothing

    # intersection point
    v1² = dot(V1,V1)
    v2² = dot(V2,V2)
    v1v2 = dot(V1,V2)
    s = ( dot(Q1-P1, V2)*v1² - dot(Q1-P1,V1)*v1v2) / (v1v2^2 - v1²*v2²)
    t = ( dot(Q1-P1, V1) + v1v2*s) / v1²

    # check if the intersection point is inside segments
    (-tol < s < 1+tol && -tol < t < 1+tol) || return nothing

    X = P1 + t*V1

    return Point(X...)
end


# Algorithm for finding loops on flat surfaces
# 1. Start with a line. If the optional list of lines is provided, the loops will only be searched inside those lines
# 2. Find all lines that share a point with the last point of the line
# 2. Set a list of visited lines (chord)
# 3. Check if the new line is coplanar with the previous lines
# 4. If the new line is coplanar, add it to the chord and go to step 2
# 6. If the chord is closed, add it to the list of loops
# 7. If the first line repeats with the in the chord, discard the chord
# 8. Continue the process after updating the list of visited lines


function findloops(line::AbstractLine; lines::Vector{<:AbstractLine}=AbstractLine[], inner=true)
    
    # get a list of lines adjacent to the last point of l
    function candidates(visited::Vector{<:AbstractLine})
        l1 = visited[end]

        if length(visited)==1
            cands = copy(line.points[end].lines)
        else
            l0 = visited[end-1]
            p1, p2 = l1.points[[1, end]]
            if p1 in l0.points
                cands = copy(p2.lines)
            else
                cands = copy(p1.lines)
            end
        end

        length(lines)>0 && intersect!(cands, lines) # warning: intersect uses hash
        filter!(!=(l1), cands)

        return cands
    end

    function findloops(visited::Vector{<:AbstractLine}, line::AbstractLine,  plane::Union{Plane,Nothing})
        if length(visited)>=2
            # todo: improve using directional graph

            idx = findfirst(==(line), visited)
            if idx==1
                visited[1].points[end] in visited[end].points && return Loop[] # closed from the wrong direction
                return [ PlaneLoop(visited) ]
            elseif idx!==nothing # discard loop (does not contain initial point)
                return PlaneLoop[]
            end

            visited[1].points[1] in line.points[[1,end]] && return [ PlaneLoop([visited; line]) ] # checking initial point
        end

        visited = [visited; line] # make a new list

        if plane===nothing && length(visited)>1
            points = [visited[1].points; line.points]
            if !colinear(points) 
                if coplanar(points)
                    plane = Plane([visited[1].points; line.points])
                else
                    return Loop[]
                end
            end
        end

        function checkplane(l::AbstractLine)
            plane===nothing && return true
            coplanar(plane, l) || return false # check if line is coplanar
            inner && return true
            # check lines with less than 2 linked surfaces:
            surfs = collect(Set(s for s in l.surfaces))
            filter!(s -> s isa PlaneFace, surfs)
            filter!(s -> s.plane==plane, surfs)
            return length(surfs)<2
        end
        
        cands = candidates(visited)
        filter!(checkplane, cands)


        # check if loop is closed
        # if length(visited)>=2
        #     line1 = visited[1]
        #     idx = findfirst(==(line1), candidates)
        #     if idx !== nothing # loop found
        #         line1.points[end] in line[end].points && return Loop
        #     end
        # end

        loops = PlaneLoop[]
        for l in cands
            append!(loops, findloops(visited, l, plane))
        end

        return loops
    end

    return findloops(AbstractLine[], line, nothing)

end


function findsurface(geo::GeoModel, l::AbstractLine)
    p1 = l.points[1]
    p2 = l.points[end]
    s1 = [s for l in p1.lines for s in l.surfaces]
    s2 = [s for l in p2.lines for s in l.surfaces]
    surfs = intersect(s1, s2)

    length(surfs)==1 && return surfs[1]
    
    for s in geo.surfaces
        s isa PlaneFace || continue
        coplanar(s.plane, l) || continue
        inside(l.points, s) && return s
    end

    return nothing
end


function addline!(geo::GeoModel, p1::Point, p2::Point; n=0, tag="")
    l = Line(p1, p2)
    ll = getline(geo, l)
    ll===nothing || return ll
    len = norm(p1.coord-p2.coord)
    
    points = Point[]
    
    # check if line intersects other lines
    for li in geo.lines
        li isa Line || continue
        p = intersection(li, l)
        p === nothing && continue
        push!(points, p)
    end

    # check if line lies over control points
    for li in geo.lines
        li isa Arc || continue
        p = li.points[2]
        insegment(p, p1, p2) && push!(points, p)
    end

    # sort points    
    push!(points, p1)
    push!(points, p2)
    points = collect(Set(points))
    sort!(points, by=p->norm(p.coord-p1.coord))

    # add missing points and update references (may update loops)
    points = [ addpoint!(geo, p) for p in points ]
    
    # add lines
    all_lines = Line[]
    npts = length(points)
    for i in 1:npts-1
        p1 = points[i]
        p2 = points[i+1]
        ni = trunc(Int, len/norm(p1.coord-p2.coord)*n)
        l = addsingleline!(geo, p1, p2, n=ni, tag=tag)
        push!(all_lines, l)

        # check if l is an endline
        (length(p1.lines)==1 || length(p2.lines)==1) && continue

        s = findsurface(geo, l)

        if s!==nothing # line divides s into two surfaces

            # find all lines inside s including the border
            lines = [ l for l in geo.lines if inside(l.points, s) ]
            lines = [ l; lines ]
            
            loops = findloops(l, lines=lines, inner=true)

            if length(loops)==0
                continue
            elseif length(loops)==1
                addplanesurface!(geo, loops[1])
            elseif length(loops)==2
                splitplanesurface!(geo, s, loops...)
            else
                error("too many loops")
            end
        else
            loops = findloops(l, inner=false)

            # add external surface
            for lo in loops
                addplanesurface!(geo, lo)
            end
        end
    end

    if length(all_lines)==1
        return all_lines[1]
    else
        return all_lines
    end
end


function addpolygon!(geo::GeoModel, points::Point...; tag="")
    n = length(points)

    lines = Line[]
    for i in 1:n
        p1 = points[i]
        p2 = points[i%n+1]
        l = addline!(geo, p1, p2, tag=tag)
        push!(lines, l)
    end

    return lines
end


function addarc!(geo::GeoModel, p1::Point, p2::Point, p3::Point; n=0, tag="")
    p1 = getpoint(geo, p1)
    p2 = getpoint(geo, p2)
    p3 = getpoint(geo, p3)

    arc = Arc(p1, p2, p3, n=n, tag=tag)
    aa = getline(geo, arc)
    aa===nothing || return aa

    # TODO: look for intersections

    # add arc
    geo._id +=1
    arc.id = geo._id
    push!(geo.lines, arc)

    push!(p1.lines, arc)
    push!(p3.lines, arc)

    loops = findloops(arc, inner=false)

    # add external surface
    for lo in loops
        addplanesurface!(geo, lo)
    end

    return arc
end


export addpath!

function addpath!(geo::GeoModel, path::Path)
    for cmd in path.cmds
        # if cmd isa LineCmd
        if cmd.key == :L
            addline!(geo, cmd.points...)
            # addline!(geo, cmd.p1, cmd.p2)
        end
    end
end

addpath!(geo::GeoModel, args...; closed=false) = addpath!(geo, Path(args...; closed=closed))



export addsubpath!
function addsubpath!(geo::GeoModel, subpath::SubPath)
    geo._id +=1
    subpath.id = geo._id
    push!(geo.subpaths, subpath)
    return subpath
end

function addsubpath!(geo::GeoModel, args...; kwargs...) 
    return addsubpath!(geo, SubPath(Path(args...); kwargs...))
end



function getpoints(lo::AbstractLoop)
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


function colinear(points::Vector{Point})
    length(points)<2 && return true
    tol = 1e-8
    
    X1 = points[1].coord
    X2 = points[2].coord
    X1X2 = X2-X1

    for p in points[3:end]
        X = p.coord
        X1X = X-X1
        N = cross(X1X2, X1X)
        if norm(N) > tol 
            return false
        end
    end
    
    return true
end


function getnormal(points::Vector{Point})
    length(points)<3 && return nothing
    tol = 1e-8
    
    X1 = points[1].coord
    X2 = points[2].coord
    X1X2 = X2-X1
    
    # find a plane looking for a non-colinear point
    local X, N, k
    for p in points[3:end]
        X = p.coord
        X1X = X-X1
        N = cross(X1X2, X1X)
        if norm(N) > tol 
            k = i+1 # index to start to check if points are in plane
            break
        end
    end
    
    # test the plane at each point
    for p in points[k:end]
        X = p.coord
        dot(X-X1, N) > tol && return nothing
    end

    return normalize(N)
end


function coplanar(points::Vector{Point})
    n = length(points)
    n<3 && return true

    tol = 1e-8
    
    # find a plane
    X1 = points[1].coord
    X2 = points[2].coord
    X1X2 = X2-X1

    # look for a non-colinear point
    local X, N, k
    for (i,p) in enumerate(points[3:n])
        X = p.coord
        X1X = X-X1
        N = cross(X1X2, X1X)
        if norm(N) > tol 
            k = i+1 # index to start to check if points are in plane
            break
        end
    end
    
    # test the plane at each point
    for p in points[k:end]
        X = p.coord
        dot(X-X1, N) > tol && return false
    end

    return true
end


function coplanar(lo::Loop)
    return coplanar(getpoints(lo))
end


# loop points have to be coplanar
function getplane(lo::Loop)
    return getplane(getpoints(lo))
end


function windingnumber(polygon::Matrix{Float64}, x::Float64, y::Float64)
    n = size(polygon, 1)
    wn = 0
    for i in 1:n
        j = i%n + 1
        xi = polygon[i,1]
        xj = polygon[j,1]
        yi = polygon[i,2]
        yj = polygon[j,2]

        # check if the ray crossed line ij
        (x<=max(xi, xj) && min(yi,yj)<=y<=max(yi,yj)) || continue

        # finds if three points are in counter clock wise sequence
        ccw = (x - xi)*(yj - yi) < (y - yi)*(xj - xi)

        k = j%n + 1
        yk = polygon[k,2] # next y in case of intersection at the vertex

        if yi<y<yj && ccw
            wn += 1
        elseif yi>y>yj && !ccw
            wn -= 1
        elseif yi<yj<yk && yj==y && ccw
            wn += 1
        elseif yi>yj>yk && yj==y && !ccw
            wn -= 1
        end
    end

    # inside if winding is not zero
    return wn
end


# todo: improve centroid calculation
function insidepolygon(testpoints::Vector{Point}, points::Vector{Point}; tol=1e-8)
    coords = [ p.coord for p in points ]
    testcoords = [ p.coord for p in testpoints ]

    # Rotating points to the xy plane.
    Z = Vec3(0,0,1)
    P = Plane(points)
    N = P.normal
    θ    = acos(dot(N, Z))
    axis = cross(N, Z)
    if norm(axis)>1e-8
        axis = normalize(axis)
        R    = Quaternion(cos(θ/2), axis[1]*sin(θ/2), axis[2]*sin(θ/2), axis[3]*sin(θ/2))
        base = sum(p.coord for p in points)/length(points)
        coords = [ base + R*(p.coord-base)*conj(R) for p in points ]
        testcoords = [ base + R*(p.coord-base)*conj(R) for p in testpoints ]
    end

    polygon   = hcat(coords...)'[:,1:2]
    # check if all testpoins are in the same level as the plane from points
    for coord in testcoords
        abs(coord[3] - coords[1][3])>1e-8 && return false
    end

    # expand or shrink polygon according to tol
    center = mean(polygon, dims=1)
    polygon = center .+ (1+tol)*(polygon .- center)

    for coord in testcoords
        x, y, _ = coord
        windingnumber(polygon, x, y)==0 && return false # outside
    end

    return true
end


function insidepolygon(pt::Point, points::Array{Point,1}; tol=1e-8)
    return insidepolygon( [pt], points, tol=tol)
end


function inside(pt::Point, loop::PlaneLoop)
    points = getpoints(loop)
    return insidepolygon(pt, points, tol=-1e-8)
end


# check if loop1 is completely inside loop2
function inside(loop1::PlaneLoop, loop2::PlaneLoop)
    # check if loops are coplanar
    loop1.plane != loop2.plane && return false

    points1 = getpoints(loop1)
    points2 = getpoints(loop2)
    length(intersect(points1, points2))>0 && return false

    !insidepolygon(points1, points2, tol=-1e-8) && return false

    return true
end

function enclosed(loop1::PlaneLoop, loop2::PlaneLoop)
    loop1.plane != loop2.plane && return false

    points1 = getpoints(loop1)
    points2 = getpoints(loop2)
    testpoints = setdiff(points1, points2)
    length(testpoints)==0 && return true  # todo ???
    return insidepolygon(testpoints, points2, tol=1e-8)
end


function inside(pts::Vector{Point}, s::PlaneFace)
    for p in pts
        !coplanar(s.plane, p) && return false
    end
    
    !insidepolygon(pts, getpoints(s.loops[1])) && return false

    for lo in s.loops[2:end]
        insidepolygon(pts, getpoints(lo), tol=-1e-8) && return false
    end

    return true
end

# function addsingleplanesurface!(geo::GeoModel, loop::Loop; tag="")

#     s1 = PlaneFace(loop, tag=tag)
#     s = getsurface(geo, s1)
#     s === nothing || return s

#     # add new loop and surface
#     loop = addsingleloop!(geo, loop)
#     geo._id +=1
#     s1.id = geo._id
#     push!(geo.surfaces, s1)

#     # update edges
#     for lo in s1.loops
#         for l in lo.lines
#             push!(l.surfaces, s1)
#         end
#     end

#     # todo: add and update volumes

#     return s1
# end


function addplanesurface!(geo::GeoModel, loop::PlaneLoop; tag="")

    surf = PlaneFace(loop, tag=tag)
    s = getsurface(geo, surf)
    s === nothing || return s

    # check if loop encloses other loops and shares any side (overlapping)
    for s in geo.surfaces
        s isa PlaneFace || continue
        surf.plane==s.plane || continue
        length(intersect(loop.lines, s.loops[1].lines))>0 && enclosed(s.loops[1], loop) && return nothing
    end

    # add new loop and surface
    loop = addsingleloop!(geo, loop)
    geo._id +=1
    surf.id = geo._id
    push!(geo.surfaces, surf)

    # check if loop is inside other surfaces (set as hole) # todo: improve for hole inside hole
    for s in geo.surfaces
        s isa PlaneFace || continue
        loop.id==s.loops[1].id && continue
        if inside(loop, s.loops[1]) 
            push!(s.loops, loop)
            for v in s.volumes
                push!(v.surfaces, surf)
                push!(surf.volumes, v)
            end
        end
    end

    # check if loop encloses other loops (set holes) # todo: improve for hole inside hole
    for s in geo.surfaces
        s isa PlaneFace || continue
        loop.id==s.loops[1].id && continue
        inside(s.loops[1], loop) && push!(surf.loops, s.loops[1])
    end

    # update edges
    for lo in surf.loops
        for l in lo.lines
            push!(l.surfaces, surf)
        end
    end

    # add volume
    surface_loops = find_surface_loops(surf)
    for loop in surface_loops
        addvolume!(geo, loop.surfaces)
    end

    return surf
end

function addsurface!(geo::GeoModel, loop::Loop; tag="")

    surf = Face(loop, tag=tag)
    s = getsurface(geo, surf)
    s === nothing || return s

    # add new loop and surface
    loop = addsingleloop!(geo, loop)
    geo._id +=1
    surf.id = geo._id
    push!(geo.surfaces, surf)

    # update edges
    for lo in surf.loops
        for l in lo.lines
            push!(l.surfaces, surf)
        end
    end

    # add volume
    surface_loops = find_surface_loops(surf)
    for loop in surface_loops
        addvolume!(geo, loop.surfaces)
    end

    return surf
end


function splitplanesurface!(geo::GeoModel, s::PlaneFace, loop1::PlaneLoop, loop2::PlaneLoop)
    # add new loop
    loop1 = addsingleloop!(geo, loop1)
    s1 = PlaneFace(loop1)
    s1.tag = s.tag

    # add surface
    geo._id +=1
    s1.id = geo._id
    push!(geo.surfaces, s1)

    # remove references to s in s.lines
    for lo in s.loops
        for l in lo.lines
            filter!(!=(s), l.surfaces)
        end
    end
    
    holeloops = s.loops[2:end]

    # update s with loop2
    s2 = s
    s2.loops[1].lines = loop2.lines

    # update holes in s1 and s2
    s1.loops = [ s1.loops[1]; [ lo for lo in holeloops if inside(lo, s1.loops[1]) ] ]
    s2.loops = [ s2.loops[1]; [ lo for lo in holeloops if inside(lo, s2.loops[1]) ] ]

    # update references to surfaces in edges
    for surf in [s1, s2]
        for lo in surf.loops
            for l in lo.lines
                push!(l.surfaces, surf)
            end
        end
    end

    # update volumes in s1
    append!(s1.volumes, s.volumes)

    # add s1 to s2.volumes
    for v in s2.volumes
        push!(v.surfaces, s1)
    end

    return s1
end


function addvolume!(geo::GeoModel, surfaces::Array{<:AbstractFace,1}; tag="")

    v = Volume(surfaces, tag=tag)
    vv = getvolume(geo, v)
    vv===nothing || return v

    geo._id +=1
    v.id = geo._id
    push!(geo.volumes, v)

    # update surfaces
    for s in v.surfaces
        push!(s.volumes, v)
    end

    return v
end


export picksurface

function picksurface(geo::GeoModel, p::Point)
    return picksurface(geo, p.coord...)
end


function picksurface(geo::GeoModel, x::Real, y::Real, z::Real=0.0)
    p = Point(x, y, z)
    for s in geo.surfaces
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


function tag!(s::AbstractFace, tag::String)
    s.tag = tag
end



function is_neighbor(surf1::AbstractFace, surf2::AbstractFace)
    for l in surf1.loops[1].lines
        for s in l.surfaces
            s==surf2 && return true
        end
    end
    return false
end


# function get_candidate_surfaces(surfs::Vector{<:AbstractFace}, border::Vector{<:AbstractLine})
#     cands = AbstractFace[]

#     # get border points
#     points = Point[]
#     for l in border
#         push!(points, l.points[1])
#         push!(points, l.points[end])
#     end
#     points = collect(Set(points))

#     # get lines shared by points
#     lines = AbstractLine[]
#     for p in points
#         for l in p.lines
#             l in lines && continue
#             push!(lines, l)
#         end
#     end

#     # get candidates
#     cands = AbstractFace[]
#     for l in lines
#         for s in l.surfaces
#             s in cands && continue
#             s in surfs && continue
#             length(s.volumes)==2 && continue # surface is between two volumes
#             any(length(l.surfaces)==1 for l in s.loops[1].lines) && continue # check for dangling edges
#             push!(cands, s)l
#         end
#     end

#     return cands
# end


# function find_surface_loops(surf::AbstractFace)
#     # surf: seed surface
#     # returns a list of surfaces that define volumes
    
#     # get neighbor surfaces (outer)
#     function get_neighbor_surfaces(surf::AbstractFace)
#         surfs = AbstractFace[]
#         for l in surf.loops[1].lines
#             for s in l.surfaces
#                 s==surf && continue
#                 s in surfs && continue
#                 length(s.volumes)==2 && continue # surface is between two volumes
#                 any(length(l.surfaces)==1 for l in s.loops[1].lines) && continue # check for dangling edges
#                 push!(surfs, s)
#             end
#         end
#         return surfs
#     end

#     # check for dangling edges
#     any(length(l.surfaces)==1 for l in surf.loops[1].lines) && return SurfaceLoop[]
#     border = surf.loops[1].lines
    
#     # get candidate surfaces
#     cands = AbstractFace[]
#     for l in border
#         for s in l.surfaces
#             s in cands && continue
#             s==surf && continue
#             length(s.volumes)==2 && continue # surface is between two volumes
#             any(length(l.surfaces)==1 for l in s.loops[1].lines) && continue # check for dangling edges
#             push!(cands, s)
#         end
#     end

#     # find crown of surfaces
#     crown = AbstractFace[]
#     for s in cands
#         crown = [ s ]
#         next = nothing
#         for ns in cands
#             if is_neighbor(s, ns) 
#                 next = ns
#                 break
#             end
#         end
#         if next==s
#             break
#         end
#         if next !== nothing
#             crown = [ crown; next ]
#         end
#     end

    
    

# end

# function find_surface_loops_recursive(surf::AbstractFace)
#     # surf: seed surface
#     # returns a list of surfaces that define volumes

#     # check for dangling edges
#     any(length(l.surfaces)==1 for l in surf.loops[1].lines) && return SurfaceLoop[]

#     function findloops(visited::Vector{<:AbstractFace}, border::Vector{<:AbstractLine}, surf::AbstractFace)
#         # surface is between two volumes
#         surf.volumes==2 && return SurfaceLoop[]

#         @show "hi"
#         @show length(border)
#         border = symdiff(border, getlines(surf))
#         @show length(border)
#         visited = [visited; surf] # make a new list
#         @show length(visited)
#         length(border)==0 && return [ SurfaceLoop(visited) ] # no border lines


#         # collect candidate surfaces
#         cands = AbstractFace[]
#         for l in border
#             for s in l.surfaces
#                 s in cands && continue
#                 s==surf && continue
#                 length(s.volumes)==2 && continue # surface is between two volumes
#                 s in visited && continue
#                 any(length(l.surfaces)==1 for l in s.loops[1].lines) && continue # check for dangling edges
#                 push!(cands, s)
#             end
#         end

#         length(cands)==0 && return SurfaceLoop[] # no candidate surfaces

#         @show length(cands)
#         for s in visited
#             @show s.id
#         end
        
        
#         loops = SurfaceLoop[]
#         for s in cands
#             # append!(loops, findloops(visited, border, s))
#             sloops = findloops(visited, border, s)
#             append!(loops, sloops)
#         end
#         @show "::::::::::::"
#         @show length(loops)
#         # error()

#         return loops
#     end

#     loops = findloops(AbstractFace[], AbstractLine[], surf)
#     # @assert length(loops)==1 # todo: check for more than one volume
#     return loops

# end


function find_surface_loops(surf::AbstractFace)
    # surf: seed surface
    # returns a list of surfaces that define volumes (currently only one volume)

    # check if surface does not add a volume
    for l in surf.loops[1].lines
        length(l.surfaces)==1 && return SurfaceLoop[]
    end

    border = surf.loops[1].lines

    # get candidate surfaces adjacent to surf
    cands = AbstractFace[]
    for l in surf.loops[1].lines
        for s in l.surfaces
            s==surf && continue
            length(s.volumes)==2 && continue # surface is between two volumes
            any(length(l.surfaces)==1 for l in s.loops[1].lines) && continue # check for dangling edges
            push!(cands, s)
        end
    end

    # reorder cands
    perm = [ length(intersect(s.loops[1].lines, border)) for s in cands ]
    cands = cands[sortperm(perm)]

    @show perm
    @show sortperm(perm)

    # find two adjacent surfaces that are also adjacent to surf (improve by getting a ring from any corner)
    surfaces = AbstractFace[ surf ]
    n = length(cands)
    for i in 1:n
        si = cands[i]
        found = false
        for j in i+1:n
            sj = cands[j]
            if is_neighbor(si, sj)
                surfaces = [ surfaces; si; sj ]
                found = true
                break
            end
        end
        found && break
    end

    # @show [ length(s.loops[1].lines) for s in surfaces ]

    if length(surfaces)==1
        return SurfaceLoop[]
    end
    
    border = symdiff(surfaces[1].loops[1].lines, surfaces[2].loops[1].lines, surfaces[3].loops[1].lines)
    visited = copy(surfaces)

    println()
    
    # grow surfaces
    while length(border)>0
        
        cands = []
        for l in border
            for s in l.surfaces
                # s in surfaces && continue
                s in visited && continue
                push!(visited, s)
                push!(cands, s)
            end
        end

        # reorder cands
        perm = [ length(intersect(s.loops[1].lines, border)) for s in cands ]
        cands = cands[sortperm(perm)]

        @show perm
        @show sortperm(perm)
        @show length(cands)
        @show length(border)

        length(cands)==0 && return SurfaceLoop[]

        # add faces
        for k in 1:1
            for s in cands
                s in surfaces && continue
                comm_edges = intersect(s.loops[1].lines, border)
                add_surf = false
                if length(comm_edges) == 1
                    if length(comm_edges[1].surfaces)==2 # edge with only two faces
                        add_surf = true
                    end
                end
                if length(comm_edges) >= 2 # face that shares at least two edges
                    add_surf = true
                end

                if add_surf
                    push!(surfaces, s)
                    border = symdiff(border, s.loops[1].lines)
                end
            end
        end

    end

    return [ SurfaceLoop(surfaces) ]

end

function find_surface_loop_orig(surf::AbstractFace)
    # surf: seed surface
    # returns a list of surfaces that define a volume

    # check if surface does not add a volume
    for l in surf.loops[1].lines
        length(l.surfaces)==1 && return SurfaceLoop[]
    end

    # find the seed edge: an edge that has the less number of surfaces
    points = Set{Point}()
    for l in surf.loops[1].lines # only in the outer loop
        push!(points, l.points[1])
        push!(points, l.points[end])
        # for p in l.points
            # push!(points, p)
        # end
    end

    points = collect(points)

    nlines = [ length(p.lines) for p in points ]
    _, idx = findmin(nlines)
    point = points[idx]
    
    if length(point.lines)==2 # point on dangling face
        return SurfaceLoop[]
    end
    
    # if length(point.lines)==0
    #     @show point.id
    #     error("point without lines")
    # end

    @assert length(point.lines)==3 # corner point | todo: check for other cases

    border_lines = getlines(surf)

    # get the edge that is not in the surf border (the one that is not coplanar
    seed_edge = setdiff(point.lines, surf.loops[1].lines)[1]

    # add two surfaces from the seed edge that share an edge with surf
    newsurfs = AbstractFace[]
    for s in seed_edge.surfaces
        s == surf && continue
        slines = getlines(s)
        if length(intersect(slines, border_lines)) == 1 # todo: >=1 ?
            push!(newsurfs, s)
        end
    end

    surfaces = [ surf; newsurfs ]

    last_n_lines = 0

    while length(border_lines)>0

        # get new border lines
        border_lines = Set{AbstractLine}()
        for s in surfaces
            for l in getlines(s)
                if l in border_lines
                    filter!(!=(l), border_lines)
                else
                    push!(border_lines, l)
                end
            end
        end

        # length(border_lines)==last_n_lines && return SurfaceLoop[] # no new border lines
        # last_n_lines = length(border_lines)

        for l in border_lines
            if length(l.surfaces)==1
                # @show "dangling EDGE"
                return SurfaceLoop[] # dangling edge
            end
        end

        # find new surfaces that share at least two edges with the border
        newsurfs = Set{AbstractFace}()
        for l in border_lines
            for s in l.surfaces
                s in surfaces && continue
                if length(intersect(getlines(s), border_lines)) >= 2
                    push!(newsurfs, s)
                end
            end
        end

        append!(surfaces, collect(newsurfs))
    end

    return [ SurfaceLoop(surfaces) ]

end