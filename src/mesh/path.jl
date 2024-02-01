export Path, addinset!

abstract type PathCmd
end


struct MoveCmd<:PathCmd
    p::Point
end

Base.length(mc::MoveCmd) = 0.0
(mc::MoveCmd)(t::Float64) = mc.p.coord

struct LineCmd<:PathCmd
    p1::Point
    p2::Point
end

Base.length(lc::LineCmd) = norm(lc.p2.coord-lc.p1.coord)
(lc::LineCmd)(t::Float64) = lc.p1.coord + t*(lc.p2.coord-lc.p1.coord) 


struct BezierCmd<:PathCmd
    p1::Point
    p2::Point
    p3::Point
    p4::Point
end

function Base.length(bc::BezierCmd) 
    P1 = bc.p1.coord
    P2 = bc.p2.coord
    P3 = bc.p3.coord
    P4 = bc.p4.coord
    deriv(t) = 3*(1-t)^2*(P2-P1) + 6*(1-t)*t*(P3-P2) + 3*t^2*(P4-P3)
    ipoints = [0.5 - √(3/5)/2, 0.5, 0.5 + √(3/5)/2]
    weights = [5/18, 8/18, 5/18]
    return sum([weights[i]*norm(deriv(ipoints[i])) for i in 1:3])
end

function (bc::BezierCmd)(t::Float64) 
    P1 = bc.p1.coord
    P2 = bc.p2.coord
    P3 = bc.p3.coord
    P4 = bc.p4.coord
    (1-t)^3*P1 + 3*(1-t)^2*t*P2 + 3*(1-t)*t^2*P3 + t^3*P4
end


mutable struct Path
    cmds::Array{PathCmd}
    closed::Bool
    _T::Array{Float64,1}
end

function Path(tokens...; closed=false)
    n  = length(tokens)
    idx = 1
    cmds = PathCmd[]
    # lastpoint = Point(0,0,0)
    # firstpoint = Point(0,0,0)
    local firstpoint, lastpoint

    while idx<=n
        token = tokens[idx]
        if token==:M
            idx+=1
            cmd = MoveCmd(tokens[idx])
            lastpoint = tokens[idx]
            if idx==2
                firstpoint = tokens[idx]
            end
        elseif token==:L
            idx+=1
            points = ( lastpoint, tokens[idx] )
            cmd = LineCmd(points...)
            lastpoint = points[end]
        elseif token==:B
            idx+=1
            points = [ lastpoint, tokens[idx]... ]
            # data = [ lastpoint'; tokens[idx] ]
            # points = [ Point(data[i,:]) for i in 1:4 ]
            cmd = BezierCmd(points...)
            lastpoint = points[end]
        end
        push!(cmds, cmd)
        
        idx+=1
    end

    if closed && firstpoint!=lastpoint
        cmd = LineCmd(lastpoint, firstpoint)
        push!(cmds, cmd)
    end

    T = cumsum( length(pc) for pc in cmds )
    T = T./T[end] # normalize


    return Path(cmds, closed, T)
end

# function Path(coords::Array, closed=false)
#     n = size(coords, 2)
#     tokens = [ :M, coords[1,:] ]
#     for i in 2:n
#         push!(tokens, :L)
#         push!(cmds, Point(coords[i,:]))
#     end

#     return Path(tokens, closed=closed)
# end

function (path::Path)(t::Float64) 
    i = searchsortedfirst(path._T, t)
    if i==1
        return path.cmds[1](0.0)
    elseif i>length(path._T)
        return path.cmds[end](1.0)
    else
        tt = (t-path._T[i-1])/(path._T[i]-path._T[i-1])
        return path.cmds[i](tt)
    end
end
