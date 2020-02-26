# This file is part of FemMesh package. See copyright license in https://github.com/NumSoftware/FemMesh

include("tools/list.jl")


export TBlock
mutable struct TBlock <: AbstractBlock
    points::Array{Point,1}
    segments::Array
    h::Float64
    sizehint::Array
    shape::ShapeType
    cellshape::ShapeType
    #ns::Array{Int,1}
    #npatch::Array{Int,1}
    tag::String
    nips::Int64
    id::Int64

    #function TBlock(coords::Array{<:Real}; h::Real=NaN, ns::Array{Int,1}=Int[], npatch::Array{Int,1}=Int[], cellshape=TRI3, tag="", id=-1, nips=0)
        #this = new()
        #points = Point[ Point(coords[i,:]) for i=1:size(coords,1) ]
        #return new(points, POLYV, cellshape, h, ns, npatch, tag, nips, id)
    #end

    function TBlock(coords::Array{<:Real}, conns::Array; h::Real=NaN, sizehint=zeros(0,0), cellshape=TRI3, tag="", id=-1, nips=0)
        this = new()
        points = Point[ Point(coords[i,:]) for i=1:size(coords,1) ]
        segments = [ [ points[i] for i in con] for con in conns  ]
        if isnan(h)
            ((minx, maxx), (miny, maxy)) = extrema(coords, dims=1)
            h = min(maxx-minx, maxy-miny)/10
        end

        return new(points, segments, h, sizehint, POLYV, cellshape, tag, nips)
    end

end

function get_h(h::Float64, sizehint::Array, X::Array{Float64,1})
    # find distances to each hint point

    n, m = size(sizehint)
    n == 0 && return h
    @assert m==3

    hmin, hmax = extrema(sizehint[:,3])
    hmin = min(h, hmin)
    hmax = max(h, hmax)

    npoints = 3 # points used for interpolation
    D = zeros(n)
    for i=1:n
        D[i] = norm(X - sizehint[i,1:2])
    end
    perm = sortperm(D)
    perm = perm[1:min(npoints,n)]
    H = sizehint[perm,:]

    # inverse distance interpolation
    sh = 0.0
    sw = 0.0
    for i=1:n
        d = norm(X-H[i,1:2])
        #@show d
        d==0 && return H[i,3]
        w = (1/d)^2
        sh += H[i,3]*w
        sw += w
    end

    #@show sh, sw

    #h = min(max(sh/sw, hmin), hmax)
    #return h
    return sh/sw
end

function split_block(bl::TBlock, msh::Mesh)
    boundary = List{Point}()

    for seg in bl.segments
        #@show seg
        #@show length(boundary)
        np = length(seg)
        if np==2
            #X1 = bl.points[seg[1],:]
            #X2 = bl.points[seg[2],:]
            #X1, X2 = seg
            X1 = [seg[1].x, seg[1].y]
            X2 = [seg[2].x, seg[2].y]
            #@show X1
            #@show X2
            len = norm(X2-X1)
            U = normalize(X2-X1)
            l = 0.0
            X = X1
            h = 0.0
            #@show len
            while l<len
                p = Point(X)
                push!(boundary, p)
                h = get_h(bl.h, bl.sizehint, X)
                l += h
                #@show h
                #@show l
                X = X1 + l*U
                norm(X2-X)<0.5*h && break
            end
        else
        end
    end

    #@show boundary
    #@show length(boundary)
    #error()

    cells, points = frontal(boundary, bl.h, bl.sizehint)
    msh.cells = cells
    msh.points = points
    return
end

function split_block2(bl::TBlock, msh::Mesh)
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

    H = zeros(n, 3)
    for i=1:n
        X1 = coords[i,:]
        X2 = coords[i+1,:]
        H[i,:] = 0.5.*(X1.+X2)
        H[i,3] = norm(X1-X2)/ns[i]
    end

    cells, points = frontal(boundary, H)
    @show length(cells)
    msh.cells = cells
    msh.points = points
end



function split_block3(bl::TBlock, msh::Mesh)
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

#=
function inv_dist_interpolation(H::Array{Float64,2}, X::Array{Float64,1})
    n, m = size(H)
    @assert m==3

    sum1 = 0.0
    sum2 = 0.0
    for i=1:n
        d = norm(X-H[i,1:2])
        #w = (1/d)^2
        w = (1/d)^1.4
        #w = (1/d)^0.5
        sum1 += H[i,3]*w
        sum2 += w
    end

    return sum1/sum2
end
=#


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

        #if abs(p.x - 0.32649766058329455) < 1e-4
        #if p.x == 0.32649766058329455
            #@show length(poly)
            #@show p1
            #@show p2
            #@show p
            #@show q
            #@show o1
            #@show o2
            #@show o3
            #@show o4
        #end

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

function frontal2(boundary::List, H::Array{Float64,2})
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
        k+=1
        #check(boundary)
        debug && println()
        debug && @show k

        #k==4 && break
        #k==25 && break
        #k==41 && break
        #k==51 && break
        #k==75 && break
        #k==120 && break
        #k==179 && break
        #k==217 && break
        #k==10000 && debug && break



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

        l = inv_dist_interpolation(H, (P1+P2)/2)
        #@show l


        l01 = √((p0.x-p1.x)^2 + (p0.y-p1.y)^2)
        l12 = √((p2.x-p1.x)^2 + (p2.y-p1.y)^2)
        l23 = √((p2.x-p3.x)^2 + (p2.y-p3.y)^2)
        #l = (l12+len)/2
        #l = (l01+l12+l23)/3
        #l = (l12+l23)/2
        #l = min(l12, l23)
        #l = min(l, len)
        #l=len
        #

        #θ1 = angle(p0, p1, p2)
        #θ2 = angle(p1, p2, p3)

        if θ>4*pi/6 || θ<pi/2
        #if θ>4*pi/6
            P = 0.5*(P1+P2) + √3/2*l*N
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
            @show l
            @show l01
            @show l12
            P = P2 + l*normalize(N1+N2)
            @show P
            @show angle(p1, p2, Point(P))*180/pi
        end

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




function frontal(boundary::List, h::Float64, H::Array{Float64,2})
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
        k+=1
        #check(boundary)
        debug && println()
        debug && @show k
        debug && @show length(boundary)

        #k==2 && break
        #k==4 && break
        #k==25 && break
        #k==41 && break
        #k==51 && break
        #k==75 && break
        #k==120 && break
        #k==179 && break
        #k==230 && break
        #k==257 && break
        #k==223 && debug && break
        #k==500 && debug && break
        #k==791 && debug && break
        #k==792 && debug && break
        #k==800 && debug && break
        #k==810 && debug && break
        #k==830 && debug && break
        #k==976 && debug && break
        #k==1096 && debug && break
        #k==10000 && debug && break



        n2 = n1.next
        n3 = n2.next

        p0 = n1.prev.data
        p1 = n1.data
        p2 = n2.data
        p3 = n3.data

        θ1 = angle(p0, p1, p2)*180/pi
        θ2 = angle(p1, p2, p3)*180/pi

        if θ1<=90 # try to close p0-p2
            if !intersects(boundary, p0, p2)
                debug && @show "CLOSE-------------p0-p2"
                cell = Cell(TRI3, [p0, p1, p2])
                push!(cells, cell)
                delete!(boundary, n1)
                n1 = n2.prev
                continue
            end
        end

        if θ2<=90 # try to close p1-p3
            if !intersects(boundary, p1, p3)
                debug && @show "CLOSE-------------p1-p3"
                cell = Cell(TRI3, [p1, p2, p3])
                push!(cells, cell)
                delete!(boundary, n2)
                #@show p1
                #@show p2
                continue
            end
        end

        @show θ1+θ2
        if 180<θ1+θ2<240 && !intersects(boundary, p0, p3)
            P0 = [p0.x, p0.y]
            P1 = [p1.x, p1.y]
            P2 = [p2.x, p2.y]
            P3 = [p3.x, p3.y]
            l02 = norm(P2-P0)
            l13 = norm(P3-P1)
            debug && @show "CLOSE-------------p0-p3"
            if l02<l13
                cell = Cell(TRI3, [p0, p1, p2])
                push!(cells, cell)
                cell = Cell(TRI3, [p0, p2, p3])
                push!(cells, cell)
            else
                cell = Cell(TRI3, [p1, p2, p3])
                push!(cells, cell)
                cell = Cell(TRI3, [p0, p1, p3])
                push!(cells, cell)
            end
            delete!(boundary, n1)
            delete!(boundary, n2)
            n1 = n3
            continue
        end

        if θ2>90 # try to add point
            P1 = [p1.x, p1.y]
            P2 = [p2.x, p2.y]
            P3 = [p3.x, p3.y]
            #l = inv_dist_interpolation(H, (P1+P2)/2)
            l = get_h(h, H, (P1+P2)/2)

            if θ2<=150 # bisector in p2
                l12 = √((p1.x-p2.x)^2 + (p1.y-p2.y)^2)
                l23 = √((p3.x-p2.x)^2 + (p3.y-p2.y)^2)
                #l = (l+l12+l23)/3
                l = 0.7*l + 0.15*l12 + 0.15*l23
                N1 = normalize!(P1-P2)
                N2 = normalize!(P3-P2)
                P = P2 + l*normalize(N1+N2)
            else # normal to p1-p2
                l12 = √((p1.x-p2.x)^2 + (p1.y-p2.y)^2)
                #l = (l+l12)/2
                l = 0.7*l+0.3*l12
                N = normalize!(P2-P1)
                N = [-N[2], N[1]]
                P = 0.5*(P1+P2) + √3/2*l*N
            end

            p = Point(P)
            pc, found = closestpoint(boundary, p, 0.9*l)
            lppc = √((p.x-pc.x)^2 + (p.y-pc.y)^2)

            #@show found

            p4 = n3.next.data
            if found && pc==p4
                if !intersects(boundary, p1, p4)
                    debug && @show "ADD1------------p1-p2-pc"
                    debug && @show "ADD1------------p2-p3-pc"
                    cell = Cell(TRI3, [p1, p2, pc])
                    push!(cells, cell)
                    cell = Cell(TRI3, [p2, p3, pc])
                    push!(cells, cell)
                    delete!(boundary, n2)
                    delete!(boundary, n3)
                    continue
                end
            elseif found && (pc==p0 || pc==p3 || lppc<0.4*l)
                debug && @show "SKIP------------------"
                n1 = n1.next
                continue
            else
                if !intersects(boundary, p2, p)
                    if !intersects(boundary, p1, p)
                        debug && @show "ADD-------------p1-p2-p"
                        cell = Cell(TRI3, [p1, p2, p])
                        push!(cells, cell)
                        push!(points, p)
                        n = Node{Point}(p)
                        insert!(boundary, n2, n)
                        if θ2<=150 && !intersects(boundary, p3, p)
                            debug && @show "ADD-------------p2-p3-p"
                            cell = Cell(TRI3, [p2, p3, p])
                            push!(cells, cell)
                            delete!(boundary, n2)
                        end
                        n1=n
                        continue
                    end
                end
            end
        end

        n1 = n1.next

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
