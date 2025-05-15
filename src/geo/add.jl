export addpoint!, addline!, addpolygon!, addarc!, addflatface!, addvolume!, addblock!


function addblock!(geo::GeoModel, block::AbstractBlock)
    push!(geo.blocks, block)

    # update ndim
    if geo.ndim != 3
        ndim = any( point.coord[3] != 0.0 for point in block.points ) ? 3 : 2
        if ndim == 2
            ndim = any( point.coord[2] != 0.0 for point in block.points ) ? 2 : 1
        end
        geo.ndim = max(ndim, geo.ndim)
    end
end


function addsinglepoint!(geo::GeoModel, pt::Point)
    p = getpoint(geo, pt)
    if p!==nothing 
        pt.size>0 && (p.size = pt.size)
        pt.tag!="" && (p.tag = pt.tag)
        return p
    end
    
    geo._id +=1
    pt.id = geo._id
    push!(geo.points, pt)

    return pt
end


function addpoint!(geo::GeoModel, pt::Point)
    p = getpoint(geo, pt)
    if p!==nothing 
        pt.size>0 && (p.size = pt.size)
        pt.tag!="" && (p.tag = pt.tag)
        return p
    end

    # add point
    geo._id +=1
    pt.id = geo._id
    push!(geo.points, pt)

    geo.ndim = pt.coord[3] != 0.0 ? 3 : max(geo.ndim, 2)

    # check if point is inside any existing line
    for edge in geo.edges
        # todo: add the case for Arcs
        if edge isa Line
            p1 = edge.points[1]
            p2 = edge.points[end]
            if in_segment(pt.coord, p1.coord, p2.coord)
                split_edge(geo, edge, pt)
            end
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

    push!(geo.edges, line)
    push!(p1.edges, line)
    push!(p2.edges, line)
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

    push!(geo.edges, arc)
    push!(p1.edges, arc)
    push!(p3.edges, arc)

    return arc
end


function addloop!(geo::GeoModel, loop::Loop)
    lo = getloop(geo, loop)
    lo===nothing || return lo

    geo._id +=1
    loop.id = geo._id
    push!(geo.loops, loop)

    return loop
end


function find_face(geo::GeoModel, l::Edge)
    p1 = l.points[1]
    p2 = l.points[end]
    s1 = [s for l in p1.edges for s in l.faces]
    s2 = [s for l in p2.edges for s in l.faces]
    surfs = intersect(s1, s2)

    length(surfs)==1 && return surfs[1]
    
    for s in geo.faces
        s.flat || continue
        coplanar(s, l.points) || continue
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
    for li in geo.edges
        li isa Line || continue
        p = intersection(li, l)
        p === nothing && continue
        push!(points, p)
    end

    # check if line lies over arc control points
    for li in geo.edges
        li isa Arc || continue
        p = li.points[2]
        in_segment(p, p1, p2) && push!(points, p)
    end

    # sort points    
    push!(points, p1)
    push!(points, p2)
    points = collect(Set(points))
    sort!(points, by=p->norm(p.coord-p1.coord))

    # add missing points and update references (may update loops)
    points = [ addpoint!(geo, p) for p in points ]
    
    # add lines
    all_edges = Line[]
    npts = length(points)
    for i in 1:npts-1
        p1 = points[i]
        p2 = points[i+1]
        ni = trunc(Int, len/norm(p1.coord-p2.coord)*n)
        l = addsingleline!(geo, p1, p2, n=ni, tag=tag)
        push!(all_edges, l)

        geo._face_detection || continue # skip face detection

        # check if l is an endline
        (length(p1.edges)==1 || length(p2.edges)==1) && continue

        # Try to add a face inside another face
        s = find_face(geo, l)
        if s!==nothing # line may divides s into two faces or make a hole

            # find all lines inside s including the border
            subset = [ l for l in geo.edges if inside(l.points, s) ] # improve
            subset = [ l; subset ]
            
            loops = find_flat_loops(Dart(l), subset, exclude_inner_edges=false)

            if length(loops)==1  # hole
                if geo._hole_filling
                    addface!(geo, loops[1], outer_face=s)
                else
                    addhole!(geo, s, loops[1])
                end
            elseif length(loops)==2  # split
                split_flat_face!(geo, s, loops...)
            elseif length(loops) > 2
                error("too many loops inside a face")
            end
        end

        # Try to add an external face
        loops = find_flat_loops(Dart(l); exclude_inner_edges=true) 

        # add external face
        for lo in loops
            lo in geo.loops && continue # loop already exists
            addface!(geo, lo)
        end
    end

    if length(all_edges)==1
        return all_edges[1]
    else
        return all_edges
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
    push!(geo.edges, arc)
    
    push!(p1.edges, arc)
    push!(p3.edges, arc)
    
    if geo._face_detection
        # loops = find_flat_loops(arc, inner=false)
        loops = find_flat_loops(Dart(arc); exclude_inner_edges=true)
        
        # add external face
        for lo in loops
            addface!(geo, lo)
        end
    end

    return arc
end



function addhole!(geo::GeoModel, face::Face, loop::Loop; tag="")
    # add new loop
    loop = addloop!(geo, loop)

    push!(face.loops, loop)
                    
    # update faces in edges
    for dart in loop.darts
        push!(dart.edge.faces, face)
    end

    # remove volumes that uses face
    for v in face.volumes
        for spin in v.spins
            filter!(!=(v), spin.face.volumes)
        end
        filter!(==(v), geo.volumes)
    end
end


function addsingleface!(geo::GeoModel, loop::Loop; tag="")
    # creates a new face with the loop
    # updates edges with the new face
    # does not trigger volume detection
    # does not chekk if it is a hole
    # does not check for inner holes
end


function addface!(geo::GeoModel, loop::Loop; outer_face=nothing, tag="")
    tol = 1e-6
    
    face = Face(loop, tag=tag)
    s = getface(geo, face)
    s === nothing || return s
    
    # check if flat face should be added: TODO: this check should performed while finding valid loops 
    if loop.flat
        normal = getnormal(loop)
        p1 = loop.darts[1].points[1]
        # check if loop encloses other loops and shares any side (overlapping)
        for s in geo.faces
            s.flat || continue
            
            # check if faces are parallel
            s_normal = getnormal(s.loops[1])
            parallel(normal, s_normal) || continue 
            
            # check if faces are coplanar
            s_p1 = s.loops[1].darts[1].points[1]
            dot(s_p1.coord - p1.coord, normal) < tol || continue # planes are at the same level
            
            # face.plane==s.plane || continue
            length(intersect(loop.darts, s.loops[1].darts))>0 && enclosed(s.loops[1], loop) && return nothing
        end
    end
    
    # add new loop and face 
    loop = addloop!(geo, loop)
    geo._id +=1
    face.id = geo._id
    push!(geo.faces, face)
    geo.quiet || println("  Added $(face.flat ? "flat" : "curved") face $(face.id)")

    # update in case of flat holes
    if face.flat
        # check if loop is inside other faces (set as hole) # todo: improve for hole inside hole
        # for f in geo.faces
        #     f.flat || continue
        #     loop.id==f.loops[1].id && continue
        #     if inside(loop, f.loops[1])
        #         push!(f.loops, loop)

        #         # update faces in edges
        #         for dart in face.loops[1].darts
        #             push!(dart.edge.faces, f)
        #         end

        #         # update volumes
        #         for v in f.volumes
        #             spin = FaceSpin(face, getnormal(f.loops[1]))
        #             push!(v.spins, spin)
        #             push!(face.volumes, v)
        #         end
        #     end
        # end

        if outer_face !== nothing # face is also a hole
            push!(outer_face.loops, loop)

            # update faces in edges
            for dart in face.loops[1].darts
                push!(dart.edge.faces, outer_face)
            end

            # update volumes
            for v in outer_face.volumes
                normal = getnormal(outer_face.loops[1])
                spin = FaceSpin(face, normal)  # todo: this spin may be wrong, see below
                push!(v.spins, spin)
                # volume_is_valid(geo, v) || (v.spins[end] = FaceSpin(face, -normal)) # todo: function volume_is_valid required
                push!(face.volumes, v)
            end
        end

        # check if loop encloses other loops (set holes) 
        # this operation might be not allowed -> error if loop encloses other loop
        # todo: improve for hole inside hole
        for f in geo.faces
            f.flat || continue
            loop.id==f.loops[1].id && continue
            if inside(f.loops[1], loop) 
                push!(face.loops, f.loops[1])
                # update faces in edges
                for dart in f.loops[1].darts
                    push!(dart.edge.faces, face)
                end
            end
        end
    end

    # update edges
    for lo in face.loops
        for dart in lo.darts
            push!(dart.edge.faces, face)
        end
    end

    # add volume
    if geo._volume_detection
        normal = getnormal(face.loops[1])
        # surface_loops = find_face_loops(FaceSpin(face, normal))
        # surface_loops = find_face_loops(FaceSpin(face, -normal))
        
        surface_loops = [ find_face_loops(FaceSpin(face, normal)); find_face_loops(FaceSpin(face, -normal)) ]
        for loop in surface_loops
            addvolume!(geo, loop.spins)
        end
    end

    return face
end


function addvolume!(geo::GeoModel, spins::Vector{FaceSpin}; tag="")

    v = Volume(spins, tag=tag)
    vv = getvolume(geo, v)
    vv===nothing || return v

    geo._id +=1
    v.id = geo._id
    push!(geo.volumes, v)

    # update faces
    for s in v.spins
        push!(s.face.volumes, v)
    end

    geo.quiet || println("  Added volume $(v.id)")

    return v
end



export addpath!

function addpath!(geo::GeoModel, path::Path)
    for cmd in path.cmds
        # if cmd isa LineCmd
        if cmd.key == :L
            addline!(geo, cmd.points...)
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