
abstract type PathCmd
end

struct MoveCmd<:PathCmd
    p::Vec3
end

Base.length(mc::MoveCmd) = 0.0
(mc::MoveCmd)(t::Float64) = mc.p

struct LineCmd<:PathCmd
    p1::Vec3
    p2::Vec3
end

Base.length(lc::LineCmd) = norm(lc.p2-lc.p1)
(lc::LineCmd)(t::Float64) = lc.p1 + t*(lc.p2-lc.p1) 


struct BezierCmd<:PathCmd
    p1::Vec3
    p2::Vec3
    p3::Vec3
    p4::Vec3
end

function Base.length(bc::BezierCmd) 
    deriv(t) = 3*(1-t)^2*(bc.p2-bc.p1) + 6*(1-t)*t*(bc.p3-bc.p2) + 3*t^2*(bc.p4-bc.p3)
    ipoints = [0.5 - √(3/5)/2, 0.5, 0.5 + √(3/5)/2]
    weights = [5/18, 8/18, 5/18]
    return sum([weights[i]*norm(deriv(ipoints[i])) for i in 1:3])
end

(bc::BezierCmd)(t::Float64) = (1-t)^3*bc.p1 + 3*(1-t)^2*t*bc.p2 + 3*(1-t)*t^2*bc.p3 + t^3*bc.p4


mutable struct Path
    cmds::Array{PathCmd}
    closed::Bool
    _T::Array{Float64,1}
end

function Path(tokens...; closed=false) 
    n  = length(tokens)
    idx = 1
    cmds = PathCmd[]
    lastpoint = Vec3(0,0,0)
    firstpoint = Vec3(0,0,0)

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
            points = ( lastpoint, Vec3(tokens[idx]) )
            cmd = LineCmd(points...)
            lastpoint = points[end]
        elseif token==:B
            idx+=1
            data = [ lastpoint'; tokens[idx] ]
            points = [ Vec3(data[i,:]) for i in 1:4 ]
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

    return Path(cmds, T, closed)
end

function Path(coords::Array, closed=false)
    n = size(coords, 2)
    tokens = [ :M, coords[1,:] ]
    for i in 2:n
        push!(tokens, :L)
        push!(cmds, Vec3(coords[i,:]))
    end

    return Path(tokens, closed=closed)
end

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



# path = Path(:M, [ 0 0 0 ], :L, [ 1 1 1 ], :B, [ 2 1 1; 2 2 1; 3 2 1] )
