# This file is part of FemMesh package. See copyright license in https://github.com/NumSoftware/FemMesh

include("tools/list.jl")


export TBlock
mutable struct TBlock <: AbstractBlock
    points::Array{Point,1}
    shape::ShapeType
    cellshape::ShapeType
    h::Float64
    ns::Array{Int,1}
    npatch::Array{Int,1}
    tag::String
    nips::Int64
    id::Int64

    function TBlock(coords::Array{<:Real}; h::Real=NaN, ns::Array{Int,1}=Int[], npatch::Array{Int,1}=Int[], cellshape=TRI3, tag="", id=-1, nips=0)
        this = new()
        points = Point[ Point(coords[i,:]) for i=1:size(coords,1) ]
        return new(points, POLYV, cellshape, h, ns, npatch, tag, nips, id)
    end

end

function split_block(bl::TBlock, msh::Mesh)
    ns = bl.ns
    n = length(ns)
    coords = getcoords(bl.points)
    np = size(coords,1)
    @assert np==n

    if length(bl.npatch)==0
        bl.npatch = repeat([1], n)
    end
    @assert sum(bl.npatch) == n

    coords = [ coords; coords[1:1,:] ]
    @show coords
    boundary = List{Point}()
    idx = 1
    for np in bl.npatch
        if np==1
            X1 = coords[idx,:]
            X2 = coords[idx+1,:]
            n1 = ns[idx]
            for i=1:n1
                X = X1 .+ (X2.-X1).*((i-1)/n1)
                p = Point(X)
                push!(boundary, p)
            end
            idx +=1
        elseif np==2
            n1 = ns[idx]
            n2 = ns[idx+1]
            C = coords[idx:idx+2, :]
            shape = LIN3
            for i=1:n1
                r = (1.0/n1)*(i-1) - 1.0
                N = bl.shape.func([r])
                X = N'*C1
                p = Point(X)
                push!(boundary, p)
            end

            for i=1:n2
                r = (1.0/n1)*(i-1)
                N = bl.shape.func([r])
                X = N'*C
                p = Point(X)
                push!(boundary, p)
            end
            idx +=2
        end
    end

    #H = zeros(n, 3)
    #for i=1:n
        #X0 = coords[i-1,:]
        #X1 = coords[i,:]
        #X2 = coords[i+1,:]
        #h1 = norm(X1-X2)/ns[i-1]
#
        #H[i,:] = 0.5.*(X1.+X2)
        #H[i,3] = norm(X1-X2)/ns[i]
    #end

    H = zeros(n, 3)
    for i=1:n
        X1 = coords[i,:]
        X2 = coords[i+1,:]
        H[i,:] = 0.5.*(X1.+X2)
        H[i,3] = norm(X1-X2)/ns[i]
    end


    #for node in boundary
        #@show node
    #end

    cells, points = frontal(boundary, H)
    @show length(cells)
    msh.cells = cells
    msh.points = points
end


function inv_dist_interpolation(H::Array{Float64,2}, X::Array{Float64,1})
    n, m = size(H)
    @assert m==3

    sum1 = 0.0
    sum2 = 0.0
    for i=1:n
        d = norm(X-H[i,1:2])
        w = (1/d)^3
        sum1 += H[i,3]*w
        sum2 += w
    end

    return sum1/sum2
end


export disc_poly
function disc_poly(points::Array{Point,1}, len::Float64)
    poly = copy(points)
    push!(poly, poly[1])

    boundary = List{Point}()
    #push!(boundary, poly[1])
    for i=1:length(points)
        p1 = poly[i]
        p2 = poly[i+1]
        L  = √( (p2.x-p1.x)^2 + (p2.y-p1.y)^2 )
        n  = floor(Int, L/len)
        ΔL = L/n

        #@show p1
        #@show p2

        push!(boundary, p1)
        for j=1:n-1
            x = p1.x + (p2.x-p1.x)/L*ΔL*j
            y = p1.y + (p2.y-p1.y)/L*ΔL*j
            #@show x, y
            p = Point(x,y)
            push!(boundary, p)
        end
    end
    return boundary
end


export isinside
"""
    isinside(poly, point)

Returns if a point is inside a polygon.
"""
function isinside(poly::List{Point}, p::Point)
    ints = 0
    n1 = poly.first
    n2 = n1.next
    p1 = n1.data
    p2 = n2.data
    #@show length(poly)
    for i=1:length(poly)
        if p.y>min(p1.y,p2.y) && p.y<=max(p1.y,p2.y)
            if p1.y != p2.y
                xi = p1.x + (p2.x-p1.x)/(p2.y-p1.y)*(p.y-p1.y)
                #@show xi
                if p.x < p1.x + (p2.x-p1.x)/(p2.y-p1.y)*(p.y-p1.y)
                    ints += 1
                end
            end
        end
        n1 = n2
        n2 = n1.next
        p1 = p2
        p2 = n2.data
    end

    #@show ints
    #@show p

    return ints%2==1
end

"""
    orientation(p1, p2, p3)

Evaluates whether a sequence of three points are disposend clockwise or counterclockwise.
Returns zero if points are colinear.
"""
function orientation(p1::Point, p2::Point, p3::Point)
    return (p2.y - p1.y) * (p3.x - p2.x) - (p2.x - p1.x) * (p3.y - p2.y)
end


function onsegment(p1::Point, p2::Point, p::Point)
    return p.x <= max(p1.x, p2.x) && p.x >= min(p1.x, p2.x) && p.y <= max(p1.y, p2.y) && p.y >= min(p1.y, p2.y)
end


export intersects
"""
    intersects(poly, p1, p2)

Returns if a segment intersects a polygon.
"""
function intersects(poly::List{Point}, p::Point, q::Point)
    n2 = poly.first
    for i=1:length(poly)
        n1 = n2
        n2 = n1.next
        p1 = n1.data
        p2 = n2.data

        max(p.x, q.x) < min(p1.x, p2.x) && continue
        min(p.x, q.x) > max(p1.x, p2.x) && continue
        max(p.y, q.y) < min(p1.y, p2.y) && continue
        min(p.y, q.y) > max(p1.y, p2.y) && continue

        o1 = orientation(p, q, p1)
        o2 = orientation(p, q, p2)
        o3 = orientation(p1, p2, p)
        o4 = orientation(p1, p2, q)

        o1*o2<0 && o3*o4<0 && return true
        #o1==0.0 && onsegment(p, q, p1) && return true
        #o2==0.0 && onsegment(p, q, p2) && return true
        #o3==0.0 && onsegment(p1, p2, p) && return true
        #o4==0.0 && onsegment(p1, p2, q) && return true
    end

    return false
end


function closestpoint(poly::List{Point}, p::Point, r::Float64)
    dmin = Inf
    n2 = poly.first
    found = false
    closest = p
    #@show length(poly)
    #check(poly)
    #@show poly.first
    #@show poly.first.next
    #@show poly.last

    for node in poly
        pp = node.data

        d = √((pp.x-p.x)^2+(pp.y-p.y)^2)
        #@show pp
        #@show d

        if d<r && d<dmin
            dmin = d
            found = true
            closest = pp
        end
    end

    return closest, found
end


function discretize_polygon(coords::Array{Float64,2}; nl::Array{Int,1}=Int[], curves::Array{Int,1}=Int[])
    points = Point[]
    n,_ = size(coords)

    if length(curves==0)
        curves = repeat([1], n)
    else
        @assert sum(curves)==n
    end

    nc = length(curves)

    ci = 0
    for i=1:n
        ni = nl[i]

    end
end



function angle(x1::Float64, y1::Float64, x2::Float64, y2::Float64)
    # vertex at (0,0)
    dotp = x1*x2 + y1*y2 # dot product
    detr = x1*y2 - y1*x2 # determinant
    return atan(detr, dotp)
end

function angle(x1::Float64, y1::Float64, x2::Float64, y2::Float64, x3::Float64, y3::Float64)
    # vertex at (x2,y2)
    dotp = (x3-x2)*(x1-x2) + (y3-y2)*(y1-y2) # dot product
    detr = (x3-x2)*(y1-y2) + (y3-y2)*(x1-x2) # determinant
    return atan(detr, dotp) # angle 321
end

function angle(p1::Point, p2::Point)
    # vertex at (0,0)
    #dotp = p1.x*p2.x + p1.y*p2.y # dot product
    #detr = p1.x*p2.y - p1.y*p2.x # determinant
    #return atan(detr, dotp)
    θ = atan(p2.y-p1.y, p2.x-p1.x)
    if θ<0; θ += 2*π end
    return θ
end

function angle(p1::Point, p2::Point, p3::Point)
    # vertex at p2
    #@show angle(p2,p3)*180/pi
    #@show angle(p2,p1)*180/pi
    θ = angle(p2,p1) - angle(p2,p3)
    if θ<0; θ += 2*π end
    return θ
    #dotp = (p3.x-p2.x)*(p1.x-p2.x) + (p3.y-p2.y)*(p1.y-p2.y) # dot product
    #detr = (p3.x-p2.x)*(p1.y-p2.y) + (p3.y-p2.y)*(p1.x-p2.x) # determinant
    #θ = atan(detr, dotp) # angle 321
    #if θ<0; θ+=2*π end
    #return θ
end

export frontal
function frontal(boundary::List, H::Array{Float64,2})
    nnodes = length(boundary)
    @assert nnodes >= 3

    cells = Cell[]
    points = [ node.data for node in boundary]

    n1 = boundary.first
    debug = false
    debug = true

    @show nnodes

try
    k=0
    while length(boundary)>2
        #check(boundary)
        debug && println()
        debug && @show k

        #k==3 && break
        #k==11 && break
        #k==41 && break
        #k==51 && break
        #k==75 && break
        #k==120 && break
        #k==207 && break
        #k==10000 && debug && break

        k+=1


        n2 = n1.next
        n3 = n2.next

        p0 = n1.prev.data
        p1 = n1.data
        p2 = n2.data
        p3 = n3.data
        θ = angle(p1, p2, p3)

        # try to add a triangle from p1 and p2
        #@show "TRYING TO ADD-------------"
        P1 = [p1.x, p1.y]
        P2 = [p2.x, p2.y]
        N = normalize!(P2-P1)
        N = [-N[2], N[1]]

        l = inv_dist_interpolation(H, P1)
        #@show l


        #l01 = √((p0.x-p1.x)^2 + (p0.y-p1.y)^2)
        l12 = √((p2.x-p1.x)^2 + (p2.y-p1.y)^2)
        l23 = √((p2.x-p3.x)^2 + (p2.y-p3.y)^2)
        #l = (l12+len)/2
        #l = (l01+l12+l23)/3
        #l = (l12+l23)/2
        #l = min(l12, l23)
        #l = min(l, len)
        #l=len
        #
        #if θ>5*pi/6 || θ<pi/2
            P = 0.5*(P1+P2) + √3/2*l*N
        #=
        else
            @show θ*180/pi
            P3 = [p3.x, p3.y]
            N1 = normalize!(P1-P2)
            N2 = normalize!(P3-P2)
            @show P1
            @show P2
            @show P3
            @show N1
            @show N2
            P = P2 + 0.9*min(l12, l23)*(N1+N2)/2
            @show P
            @show angle(p1, p2, Point(P))*180/pi
            #error()
        end
        =#

        p = Point(P)

        pc, found = closestpoint(boundary, p, 0.9*l)
        @show pc

        if !found && isinside(boundary, p)
            debug && @show "ADDING-------------"
            cell = Cell(TRI3, [p1, p2, p])
            push!(cells, cell)
            push!(points, p)
            n = Node{Point}(p)
            insert!(boundary, n2, n)
            n1=n
        elseif found && pc==p0
            debug && @show "CLOSE-------------p0"
            cell = Cell(TRI3, [p0, p1, p2])
            push!(cells, cell)
            delete!(boundary, n1)
            θ = angle(p0, p2, p3)
            if θ<π/2
                n1 = n2.prev
            else
                n1 = n2
            end
        elseif found && pc==p3
            debug && @show "CLOSE-------------p3"
            cell = Cell(TRI3, [p1, p2, p3])
            push!(cells, cell)
            delete!(boundary, n2)
            θ = angle(p0, p1, p3)
            if θ<π/2
                n1 = n1.prev
            end
        elseif θ<π/6 # Warnign: there is no intersection check
            debug && @show "CLOSE-------------θ"
            cell = Cell(TRI3, [p1, p2, p3])
            push!(cells, cell)
            delete!(boundary, n2)
            θ = angle(p0, p1, p3)
            if θ<π/2
                n1 = n1.prev
            end
        else
            debug && @show "SKIPPING-------------"
            n1 = n1.next
        end

        debug && @show length(boundary)
    end
    debug && println("Triangulation finished!!!")
    debug && @show length(cells)
    debug && @show length(points)

catch e
    printstyled(e, color=:red, bold=true)
    println()
end

    return cells, points
end
