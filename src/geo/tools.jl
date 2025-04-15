

function flip(dart::Dart)
    return Dart(dart.edge, forward=!dart.forward)
end


function dir(dart::Dart)
    X2 = dart.points[end].coord
    X1 = dart.points[1].coord
    return normalize(X2-X1)
end



function getedges(s::Face)
    edges = Edge[]
    for loop in s.loops
        for dart in loop.darts
            push!(edges, dart.edge)
        end
    end
    return edges
end


function getpoints(lo::Loop)
    points = Point[]

    for dart in lo.darts
        push!(points, dart.points[1])
    end

    return points

end


function getpoints(s::Face)
    points = Set{Point}()
    for edge in getedges(s)
        for point in edge.points
            push!(points, point)
        end
    end
    return collect(points)
end



function getpoints(v::Volume)
    points = Set{Point}()
    for spin in v.spins
        for point in getpoints(spin.face)
            push!(points, point)
        end
    end
    return collect(points)
end



function is_neighbor(surf1::Face, surf2::Face)
    for dart in surf1.loops[1].darts
        for s in dart.edge.faces
            s==surf2 && return true
        end
    end
    return false
end




function in_segment(X, X1, X2)
    tol = 1e-6
    X1X = X-X1
    X1X2 = X2-X1
    norm(cross(X1X2, X1X)) > tol && return false # check colinearity
    dot1 = dot(X1X2, X1X)
    dot2 = dot(X1X2, X1X2)
    return tol < dot1 < dot2-tol # check if inside segment within tol
end

function in_segment(p::Point, p1::Point, p2::Point)
    return in_segment(p.coord, p1.coord, p2.coord)
end


@enum(GeoFindStatus,
    OUTSIDE  = 0,
    INSIDE   = 1,
    ATPOINT  = 2,
    ATBORDER = 3,
)

function point_in_edge(p::Point, l::Edge)

    p in l.points[[1,end]] && return ATPOINT

    tol = 1e-6
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


function coplanar(face::Face, points::Vector{Point})
    face.flat || return false
    tol = 1e-6 # reduced tolerance 

    normal = getnormal(face.loops[1])
    p1 = face.loops[1].darts[1].points[1]

    for p in points # check if every point is in the plane
        abs(dot(p.coord - p1.coord, normal)) > tol && return false
    end

    return true
end


function split_edge(geo::GeoModel, edge::Edge, point::Point)

    p1 = edge.points[1]
    p2 = edge.points[end]

    if edge isa Line
        n = div(edge.n,2)
        edge1 = addsingleline!(geo, p1, point, n=n, tag=edge.tag)
        edge2 = addsingleline!(geo, point, p2, n=n, tag=edge.tag)
    elseif edge isa Arc
        error("split_edge > Arc : Not implemented")
    end

    newedges = Edge[ edge1, edge2 ]
    
    # update point edges
    point.edges = newedges

    # remove original edge
    filter!(!=(edge), geo.edges) 

    #remove references to edge in points
    filter!(!=(edge), p1.edges)
    filter!(!=(edge), p2.edges)

    # update face references
    edge1.faces = copy(edge.faces)
    edge2.faces = copy(edge.faces)

    # update loops that contained the original edge
    loops = [ loop for face in edge.faces for loop in face.loops ]
    for loop in loops
        edges = [ dart.edge for dart in loop.darts ]
        idx = findfirst(==(edge), edges)
        idx===nothing && continue # skip if not found (e.g. if loop is a hole)
        dart = loop.darts[idx]
        
        if dart.forward
            new_darts = [ Dart(edge1, forward=true), Dart(edge2, forward=true) ]
        else
            new_darts = [ Dart(edge2, forward=false), Dart(edge1, forward=false) ]
        end

        loop.darts = [ loop.darts[1:idx-1]; new_darts; loop.darts[idx+1:end] ]
    end

    return newedges
end





function intersection(l1::Edge, l2::Edge)
    # Assumes:
    # L1: P = P1 +V1*t
    # L2: Q = Y2 +V2*s
    # V1 = P2-P1
    # V2 = Q2-Q1

    tol = 1e-6
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


# check if two vectors are parallel
function parallel(V1::AbstractArray{Float64}, V2::AbstractArray{Float64})
    tol = 1e-6
    return norm(cross(V1,V2)) < tol
end


# check if two darts are parallel
function parallel(dart1::Dart, dart2::Dart)
    V1 = dir(dart1)
    V2 = dir(dart2)
    return parallel(V1, V2)
end





function colinear(points::Vector{Point})
    length(points)<2 && return true
    tol = 1e-6
    
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


function getnormal(loop::Loop)
    C = getcoords(loop.points)
    
    if !loop.flat # compute the projection to the least squares plane
        # Get centered coordinates
        center = vec(mean(C, dims=1))
        C = C .- center'
        
        # Perform SVD on the centered points
        _, _, Vt = svd(C)
        n = Vt[:, end] # least-squares plane normal (unit)
        
        # Project points onto the plane
        C = C .- (C * n) * n'
    end
        
    # Apply Newell’s method to projected (coplanar) points
    nx = 0.0
    ny = 0.0
    nz = 0.0
    n = size(C, 1)
    for i in 1:n
        X1 = C[i, :]
        X2 = C[mod1(i+1, n), :]  # wrap around

        nx += (X1[2] - X2[2]) * (X1[3] + X2[3])
        ny += (X1[3] - X2[3]) * (X1[1] + X2[1])
        nz += (X1[1] - X2[1]) * (X1[2] + X2[2])
    end

    return normalize(Vec3(nx, ny, nz))

end


function coplanar(points::Vector{Point})
    n = length(points)
    n<=3 && return true

    tol = 1e-6
    
    # find a plane
    X1 = points[1].coord
    X2 = points[2].coord
    X1X2 = X2-X1

    # look for a non-colinear point
    local X, N
    k = 1
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
        abs(dot(X-X1, N)) > tol && return false
    end

    return true
end

# Checks if a set of points are coplanar to the plane defined by a point and a normal
function coplanar(point::Point, normal::AbstractArray{Float64}, points::Vector{Point})
    tol = 1e-6 # reduced tolerance 

    for p in points # check if every point is in the plane
        abs(dot(p.coord - point.coord, normal)) > tol && return false
    end

    return true
end


# check if darts are coplanar (include arcs)
function coplanar(dart1::Dart, dart2::Dart)

    points = [ dart1.edge.points; dart2.edge.points ]
    # points = [ dart1.edge.points[[1,end]]; dart2.edge.points[[1,end]] ]
    # dart1.edge isa Arc && push!(points, dart1.edge.extrapoints[2])
    # dart2.edge isa Arc && push!(points, dart2.edge.extrapoints[2])
    return coplanar(points)

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

    # coplanar(points) || return false

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


function inside(pt::Point, loop::Loop)
    points = getpoints(loop)
    return insidepolygon(pt, points, tol=-1e-8)
end


# check if loop1 is completely inside loop2
function inside(loop1::Loop, loop2::Loop)
    @assert loop1.flat==loop2.flat==true
    tol = 1e-6

    # check if loops are parallel
    n1 = getnormal(loop1)
    n2 = getnormal(loop2)
    parallel(n1, n2) || return false

    # check if loops are coplanar
    p1 = loop1.darts[1].points[1]
    p2 = loop2.darts[1].points[1]
    dot(p1.coord - p2.coord, n1) < tol || return false

    points1 = getpoints(loop1)
    points2 = getpoints(loop2)
    length(intersect(points1, points2))>0 && return false

    !insidepolygon(points1, points2, tol=-1e-8) && return false

    return true
end


function enclosed(loop1::Loop, loop2::Loop)
    @assert loop1.flat==loop2.flat==true

    tol = 1e-6

    # check if loops are parallel
    n1 = getnormal(loop1)
    n2 = getnormal(loop2)
    parallel(n1, n2) || return false

    # check if loops are coplanar
    p1 = loop1.darts[1].points[1]
    p2 = loop2.darts[1].points[1]
    dot(p1.coord - p2.coord, n1) < tol || return false
    # loop1.plane != loop2.plane && return false

    points1 = getpoints(loop1)
    points2 = getpoints(loop2)
    testpoints = setdiff(points1, points2)
    length(testpoints)==0 && return true  # todo ???
    return insidepolygon(testpoints, points2, tol=1e-8)
end


function inside(points::Vector{Point}, s::Face)
    @assert s.flat
    tol = 1e-6
    normal = getnormal(s.loops[1])
    p1 = s.loops[1].darts[1].points[1]
    
    # check if all points are coplanar with s
    coplanar(s, points) || return false
    # for p in points
    #     dot(p.coord - p1.coord, normal) > tol && return false
    #     #x !coplanar(s.plane, p) && return false
    # end
    
    !insidepolygon(points, getpoints(s.loops[1])) && return false

    for lo in s.loops[2:end]
        insidepolygon(points, getpoints(lo), tol=-1e-8) && return false
    end

    return true
end