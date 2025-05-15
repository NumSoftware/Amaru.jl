
function revolve!(
    geo::GeoModel, 
    line::Line; 
    base::Vector{<:Real} = Float64[],
    axis::Vector{<:Real} = Float64[],
    angle::Real = NaN,
)
    @check length(base)==3
    @check length(axis)==3

    axis = Vec3(normalize(axis))
    base = Vec3(base)

    θ = angle*pi/180
    R = Quaternion(cos(θ/2), axis[1]*sin(θ/2), axis[2]*sin(θ/2), axis[3]*sin(θ/2))

    p1 = line.points[1]
    p2 = line.points[2]
    
    C3 = base + R*(p1.coord - base)*conj(R)
    C4 = base + R*(p2.coord - base)*conj(R)
    p3 = Point(C3, size=p1.size)
    p4 = Point(C4, size=p2.size)

    p1==p3 && p2==p4 && return nothing, line # line is on the axis
    
    # center points
    C5 = base + ( dot(p1.coord - base, axis) )*axis
    C6 = base + ( dot(p2.coord - base, axis) )*axis
    p5 = addsinglepoint!(geo, Point(C5))
    p6 = addsinglepoint!(geo, Point(C6))

    line2 = addline!(geo, p3, p4)
    # line2 = addsingleline!(geo, p3, p4)

    edges = Edge[ ]
    if p1!=p3
        arc = addarc!(geo, p1, p5, p3)
        push!(edges, arc)
    end

    push!(edges, line)

    if p2!=p4
        arc = addarc!(geo, p2, p6, p4)
        push!(edges, arc)
    end

    push!(edges, line2)

    flat = p5==p6 # same center -> flat face
    loop = Loop(edges..., flat=flat)
    s = addface!(geo, loop)
    
    return s, line2
end


function revolve!(
    geo::GeoModel, 
    face::Face; 
    base::Vector{<:Real} = Float64[],
    axis::Vector{<:Real} = Float64[],
    angle::Real = NaN,
)
    @check length(base)==3
    @check length(axis)==3
    @check -360 <= angle <= 360

    steps = ceil(Int, abs(angle)/120)
    angle_inc = angle/steps
    tag = face.tag

    # disable automatic volume detection
    geo._volume_detection = false
    # disable automatic face detection
    geo._face_detection = false
    # disable hole filling
    geo._hole_filling = false

    # revolve all lines
    final_lines = Edge[]
    new_faces = Face[]
    for lo in face.loops
        for edge in lo.edges
            line = edge
            for i in 1:steps
                s, line = revolve!(geo, edge, base=base, axis=axis, angle=angle_inc)
                s === nothing && break # line is on the axis
                edge = line # get the opposite line
                push!(new_faces, s)
            end
            push!(final_lines, line)
        end
    end

    # find most distant point to the axis (to be used later)
    dist = 0.0
    idx = 0
    points = getpoints(face)
    for (i,p) in enumerate(points)
        V = p.coord - base
        d1 = dot(V, axis)
        d2 = dot(V, V) - d1*d1
        if d2>dist
            dist = d2
            idx = i
        end
    end
    dpoint = points[idx]
    
    # face normal
    normal = getnormal(face.loops[1])

    if mod(angle,360) == 0
        # remove face from geo
        delete!(geo, face)

        # find most distant dart in face
        dist = 0.0
        idx = 0
        d_dart = nothing
        for (i,dart) in enumerate(face.loops[1].darts)
            edge = dart.edge
            X1 = edge.points[1].coord
            X2 = edge.points[end].coord
            V = 0.5*(X1 + X2) - base
            d1 = dot(V, axis)
            d2 = dot(V, V) - d1*d1
            if d2>dist && dpoint in edge.points
                dist = d2
                idx = i
                d_dart = dart
            end
        end

        # normal of the original face
        normal = getnormal(face.loops[1])

        # update face with one that contains d_dart.edge
        for f in new_faces
            for d in f.loops[1].darts
                if d.edge == d_dart.edge
                    face = f
                    break
                end
            end
        end

        X1 = d_dart.points[1].coord
        X2 = d_dart.points[end].coord
        V  = X2 - X1
        Nr = cross(V, normal) # points outwards
        
        normal = getnormal(face.loops[1]) # normal of the new seed face
        dot(normal, Nr) < 0.0 && (normal = -normal)
    else
        # add lid face
        loops = find_flat_loops(Dart(final_lines[1]), final_lines, exclude_inner_edges=true)
        loop = addloop!(geo, loops[1])
        addface!(geo, loop)

        # find normal
        V = dpoint.coord - base
        Nr = cross(axis, V)
        normal = getnormal(face.loops[1])
        dot(normal, Nr) > 0.0 && (normal = -normal)
    end

    # add volume
    # loop = find_face_loops(FaceSpin(face, normal))[1]
    # addvolume!(geo, loop.spins)
    loops = find_face_loops(FaceSpin(face, normal))
    if Base.length(loops)>0
        addvolume!(geo, loops[1].spins, tag=tag)
    end

    # restore automatic volume detection
    geo._volume_detection = true
    # restore automatic face detection
    geo._face_detection = true
    # restore hole filling
    geo._hole_filling = true
end


function revolve!(
    geo::GeoModel, 
    faces::Vector{Face};
    base::Vector{<:Real} = Float64[],
    axis::Vector{<:Real} = Float64[],
    angle::Real = NaN,
)
    @check length(base)==3
    @check length(axis)==3
    @check -360 <= angle <= 360

    # revolve all lines
    faces = copy(faces) # detach from geo.surfaces since the latter gets updated
    for s in faces
        revolve!(geo, s, base=base, axis=axis, angle=angle)
    end
end