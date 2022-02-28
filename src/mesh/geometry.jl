export Point, Line, Polygon, PolygonMesh

mutable struct Point<:AbstractPoint
    id::Int
    coord::Vec3
    tag::String
    function Point(x::Real, y::Real=0.0, z::Real=0.0; tag="", id=-1)
        x = round(x, digits=8)
        y = round(y, digits=8)
        z = round(z, digits=8)
        return new(id, Vec3(x,y,z), "")
    end

    function Point(X::AbstractArray{<:Real}; tag="", id=-1)
        @assert length(X) in (1,2,3)
        return Point(X..., tag=tag, id=id)
    end
end

Base.hash(p::Point) = hash( (p.coord.x, p.coord.y, p.coord.z) )


# Index operator for an collection of points
function Base.getindex(points::Array{Point,1}, s::String)
    R = [ point for point in points if point.tag==s ]
    sort!(R, by=point->sum(point.coord))
end

# Index operator for an collection of points
function Base.getindex(points::Array{Point,1}, filter_ex::Expr)
    R = Point[]
    for point in points
        x, y, z = point.coord
        eval_arith_expr(filter_ex, x=x, y=y, z=z) && push!(R, point)
    end

    sort!(R, by=point->sum(point.coord))
    return R
end



mutable struct Line
    points::Array{Point,1}
end
Base.hash(l::Line) = sum(hash(p) for p in l.points)

mutable struct Polygon
    points::Array{Point,1}
    lines::Array{Line,1}

    function Polygon(points::Array{Point,1})
        npoints = length(points)
        lines = Line[]
        for (p1,p2) in zip(points, points[[2:npoints;1]])
            push!(lines, Line([p1,p2]))
        end
        return new(points, lines)
    end

    function Polygon(coords::Array{<:Real,2})
        points = Point[ Point(row) for row in eachrow(coords) ]
        return Polygon(points)
    end
end
Base.hash(poly::Polygon) = sum(hash(p) for p in poly.points)

# Index operator for a collection of polygons using an expression
function Base.getindex(polys::Array{Polygon,1}, filter_ex::Expr)
    length(polys)==0 && return Element[]

    #points = polys[:points] # Todo
    points = collect(Set( p for poly in polys for p in poly.points ))
    rdict = Dict{Point,Bool}()
    for (i,point) in enumerate(points)
        x, y, z = point.coord
        rdict[point] = eval_arith_expr(filter_ex, x=x, y=y, z=z)
    end

    R = Polygon[]
    for poly in polys
        all( rdict[point] for point in poly.points ) && push!(R, poly)
    end
    return R
end

# function to check if two segments intersect
function intersect(l1::Line, l2::Line)

end


function Base.contains(poly::Polygon, point::Point)

    # Get avg plane

    # coordinates
    eps = 1e-8
    npoints = length(poly.points)
    C = [  p.coord[j] for p in poly.points, j in 1:3 ]

    # display(C)

    C .+= eps
    I = ones(npoints)
    N = pinv(C)*I # best fit normal
    N = normalize(N)

    # project to plane
    v1 = normalize(poly.points[2].coord - poly.points[1].coord)
    v2 = normalize(cross(N, v1))
    v1 = cross(v2, N)

    R = [v1'; v2'; N']
    XY = (C*R')[:,1:2]
    XY .-= eps
    XY = round.(XY, digits=8)
    # display(XY)

    X = R*point.coord
    x, y = X[1], X[2]

    # find if point is inside
    ints = 0

    for i in 1:npoints
        x1, y1 = XY[i,1], XY[i,2]

        # check if point is equal to vertex
        abs(x1-x)<eps && abs(y1-y)<eps && return true

        # get second point of segment
        if i!=size(XY,1)
            x2, y2 = XY[i+1,1], XY[i+1,2]    
        else
            x2, y2 = XY[1,1], XY[1,2]    
        end

        y1==y2 && continue

        if y>=min(y1,y2) && y<=max(y1,y2)
            # check if point is contained in line
            abs((x2-x1)/(y2-y1) - (x2-x)/(y2-y)) < eps && return true

            xi = x1 + (x2-x1)/(y2-y1)*(y-y1)

            if xi > x 
                ints += 1
            end
        end

    end

    return ints%2==1

end


mutable struct PolygonMesh
    points::Array{Point,1}
    lines::Array{Line,1}
    polygons::Array{Polygon,1}

    function PolygonMesh(polygons::Array{Polygon,1})
        pdict = Dict{UInt,Point}( hash(p)=>p for poly in polygons for p in poly.points )
        points = collect(values(pdict))
        for poly in polygons
            poly.points = Point[ pdict[hash(p)] for p in poly.points ]
        end

        ldict = Dict{UInt,Line}( hash(l)=>l for poly in polygons for l in poly.lines )
        lines = collect(values(ldict))
        for line in lines
            line.points = Point[ pdict[hash(p)] for p in line.points ]
        end

        for poly in polygons
            poly.lines = Line[ ldict[hash(l)] for l in poly.lines ]
            #poly.points = Point[ pdict[hash(p)] for p in poly.points ]
            #for line in poly.lines
                #line.points = Point[ pdict[hash(p)] for p in line.points ]
            #end
        end

        return new(points, lines, polygons)
    end
end



# Tag functions
for T in (Point, Line, Polygon)
    @eval begin
        tag!(object::$T, tag::String) = (object.tag = tag)
        tag!(objects::Array{<:$T,1}, tag::String) = for object in objects; object.tag = tag end
    end
end




function extrude(line::Line; dir=[0,0,1], length=1.0)
    p1, p2 = line.points
    p3 = Point(p1.coord .+ length.*dir)
    p4 = Point(p2.coord .+ length.*dir)
    return Polygon([p1,p2,p4,p3])
end

function extrude(poly::Polygon; dir=[0,0,1], length=1.0)
    polygons = [ poly ]

    # lateral polygons
    for line in poly.lines
        push!(polygons, extrude(line, dir=dir, length=length))
    end

    # back polygon
    backpoints = [ Point(p.coord .+ length.*dir) for p in poly.points ] 
    push!(polygons, Polygon(backpoints))

    return PolygonMesh(polygons)

end

function getcoords(polym::PolygonMesh)
    return Float64[ p.coord[j] for p in polym.points, j=1:3 ]
end

function getcoords(points::Array{Point,1})
    return Float64[ p.coord[j] for p in points, j=1:3 ]
end
