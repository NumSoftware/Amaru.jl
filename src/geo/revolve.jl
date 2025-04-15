
function revolve!(
    geo::GeoModel, 
    line::Edge; 
    base::Vector{<:Real} = Float64[],
    axis::Vector{<:Real} = Float64[],
    angle::Real = NaN,
)
    @check length(base)==3
    @check length(axis)==3
    @check -180 < angle < 180

    axis = Vec3(normalize(axis))
    base = Vec3(base)

    θ = angle*pi/180
    R = Quaternion(cos(θ/2), axis[1]*sin(θ/2), axis[2]*sin(θ/2), axis[3]*sin(θ/2))
    
    C1 = line.points[1].coord
    C2 = line.points[2].coord
    C3 = base + R*(C1 - base)*conj(R)
    C4 = base + R*(C2 - base)*conj(R)

    sz1 = line.points[1].size
    sz2 = line.points[2].size

    p1 = Point(C1, size=sz1)
    p2 = Point(C2, size=sz2)
    p3 = Point(C3, size=sz1)
    p4 = Point(C4, size=sz2)

    if p1==p3 && p2==p4 # line is on the axis
        return nothing
    end
    
    # center points
    C5 = base + ( dot(C1-base, axis) )*axis
    C6 = base + ( dot(C2-base, axis) )*axis
    p5 = addsinglepoint!(geo, Point(C5)) 
    p6 = addsinglepoint!(geo, Point(C6))

    l2 = addline!(geo, p3, p4)
    # l2 = addsingleline!(geo, p3, p4)

    if p1==p3 # p1 is on the axis
        a2 = addarc!(geo, p2, p6, p4)
        if p5==p6 # same center
            # flat face
            loop = Loop(line, a2, l2, flat=true)
            s = addface!(geo, loop)
        else
            loop = Loop(line, a2, l2, flat=false)
            s = addface!(geo, loop)
        end    
    elseif p2==p4 # p2 is on the axis
        a1 = addarc!(geo, p1, p5, p3)
        if p5==p6 # same center
            # flat face
            loop = Loop(line, a1, l2, flat=true)
            s = addface!(geo, loop)
        else
            loop = Loop(line, a1, l2, flat=false)
            s = addface!(geo, loop)
        end
    else # p1 and p2 are not on the axis
        a1 = addarc!(geo, p1, p5, p3)
        a2 = addarc!(geo, p2, p6, p4)
        if p5==p6 # same center
            loop = Loop(line, a1, l2, a2, flat=true)
            s = addface!(geo, loop)
        else
            loop = Loop(line, a1, l2, a2, flat=false)
            s = addface!(geo, loop)
        end
    end

    return s
end


function revolve!(
    geo::GeoModel, 
    face::Face; 
    base::Vector{<:Real} = Float64[],
    axis::Vector{<:Real} = Float64[],
    angle   ::Real = NaN,
)
    @check length(base)==3
    @check length(axis)==3
    @check -180 < angle < 180

    geo._volume_detection = false 

    # revolve all lines
    for lo in face.loops
        for edge in lo.edges
            revolve!(geo, edge, base=base, axis=axis, angle=angle)
        end
    end

    # get most distant point to the axis
    dist = 0.0
    idx = 0
    points = getpoints(face)
    for (i,p) in enumerate(points)
        V = p.coord - base
        d1 = dot(V, axis)
        d = (dot(V, V) - d1*d1)^0.5
        if d>dist
            dist = d
            idx = i
        end
    end
    V = points[idx].coord - base
    Nr = cross(axis, V)

    # find normal
    normal = getnormal(face.loops[1])
    dot(normal, Nr) > 0.0 && (normal = -normal)

    # loop = find_face_loops(FaceSpin(face, normal))[1]
    # addvolume!(geo, loop.spins)
    loops = find_face_loops(FaceSpin(face, normal))
    if Base.length(loops)>0
        addvolume!(geo, loops[1].spins)
    end

    geo._volume_detection = true 
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
    @check -180 < angle < 180

    # revolve all lines
    faces = copy(faces) # detach from geo.surfaces since the latter gets updated
    for s in faces
        revolve!(geo, s, base=base, axis=axis, angle=angle)
    end
end