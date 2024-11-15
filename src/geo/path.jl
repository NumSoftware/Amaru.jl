export Path, addinset!

abstract type PathCmd
end


struct MoveCmd<:PathCmd
    p::Point
end

Base.length(mc::MoveCmd) = 0.0
startpoint(mc::MoveCmd) = mc.p
endpoint(mc::MoveCmd) = mc.p
(mc::MoveCmd)(t::Float64) = mc.p.coord

struct LineCmd<:PathCmd
    p1::Point
    p2::Point
end

Base.length(lc::LineCmd) = norm(lc.p2.coord-lc.p1.coord)
startpoint(lc::LineCmd) = lc.p1
endpoint(lc::LineCmd) = lc.p2
(lc::LineCmd)(t::Float64) = lc.p1.coord + t*(lc.p2.coord-lc.p1.coord) 


struct ArcCmd<:PathCmd
    p1::Point
    p2::Point # center
    p3::Point
end

startpoint(ac::ArcCmd) = ac.p1
endpoint(ac::ArcCmd) = ac.p3

function Base.length(ac::ArcCmd)
    X1 = ac.p1.coord
    X2 = ac.p2.coord
    X3 = ac.p3.coord
    r = norm(X2-X1)
    θ = acos(clamp(dot(X1-X2, X3-X2)/r^2, -1, 1))
    len = r*θ
    
    return len
end

function (bc::ArcCmd)(t::Float64)
    X1 = ac.p1.coord
    X2 = ac.p2.coord
    X3 = ac.p3.coord
    axis = cross(X1-X2, X3-X2)
    axis = axis/norm(axis)
    r = norm(X2-X1)
    θ = t*acos(clamp(dot(X1-X2, X3-X2)/r^2, -1, 1))
    R = Quaternion(cos(θ/2), axis[1]*sin(θ/2), axis[2]*sin(θ/2), axis[3]*sin(θ/2))
    P = X1 + R*(X2-X1)*conj(R)
    return P
end

struct BezierCmd<:PathCmd
    p1::Point
    p2::Point
    p3::Point
    p4::Point
end

startpoint(bc::BezierCmd) = bc.p1
endpoint(bc::BezierCmd) = bc.p4

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


function path_from_numbers(tokens...; closed=false)

    n  = length(tokens)
    idx = 1
    cmds = PathCmd[]
    local startpoint, endpoint
    @show tokens

    kidx = 1
    while kidx !== nothing
        key = tokens[kidx]
        nidx = kidx+1
        kidx = findfirst( i->isa(i, Symbol), tokens[idx+1:end] )
        if kidx === nothing
            data = tokens[nidx:end]
        else
            data = tokens[nidx:kidx-1]
        end

        if key==:M
            p = Point(data...)
            cmd = MoveCmd(p)
            endpoint = p
            if idx==2
                startpoint = p
            end
        elseif key==:L
            p2 = Point(data...)
            cmd = LineCmd(endpoint, p2)
            endpoint = p2
        elseif key==:B
            if length(data)==6
                p2 = Point(data[1:2]...)
                p3 = Point(data[3:4]...)
                p4 = Point(data[5:6]...)
            elseif length(data)==9
                p2 = Point(data[1:3]...)
                p3 = Point(data[4:6]...)
                p4 = Point(data[7:9]...)
            else
                error("Invalid number of points for Bezier curve")
            end
            cmd = BezierCmd(endpoint, p2, p3, p4)
            endpoint = p4
        end
        push!(cmds, cmd)
    end

    if closed && firstpoint!=endpoint
        cmd = LineCmd(endpoint, startpoint)
        push!(cmds, cmd)
    end

    if !closed && startpoint===endpoint
        closed = true
    end

    T = cumsum( length(pc) for pc in cmds )
    T = T./T[end] # normalize


    return Path(cmds, closed, T)
end

function Path(tokens...; closed=false)
    any( t->isa(t, Point), tokens ) || return path_from_numbers(tokens...; closed=closed)

    n  = length(tokens)
    idx = 1
    cmds = PathCmd[]
    local startpoint, endpoint

    while idx<=n
        token = tokens[idx]
        if token==:M
            cmd = MoveCmd(tokens[idx+1])
            endpoint = cmd.p
            if idx==1
                startpoint = cmd.p
            end
            idx+=2
        elseif token==:L
            cmd = LineCmd(endpoint, tokens[idx+1])
            endpoint = cmd.p2
            idx+=2
        elseif token==:A
            cmd = ArcCmd(endpoint, tokens[idx+1], tokens[idx+2])
            endpoint = cmd.p3
            idx+=3
        elseif token==:C
            cmd = BezierCmd(endpoint, tokens[idx+1:idx+3]...)
            endpoint = cmd.p4
            idx+=4
        else
            error("Invalid token $token")
        end
        push!(cmds, cmd)
        
    end

    if closed && firstpoint!=endpoint
        cmd = LineCmd(endpoint, startpoint)
        push!(cmds, cmd)
    end

    if !closed && startpoint===endpoint
        closed = true
    end

    T = cumsum( length(pc) for pc in cmds )
    T = T./T[end] # normalize


    return Path(cmds, closed, T)
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


mutable struct SubPath<:GeoEntity
    id::Int
    path::Path
    closed ::Bool
    embedded ::Bool
    shape::CellShape
    tag      ::String
    jointtag ::String
    tipjointtag::String
    tipjoint::Symbol

    function SubPath(
        path::Path;
        closed::Bool  = false,
        embedded::Bool  = false,
        shape   ::CellShape = LIN3,
        tag     ::String  = "",
        jointtag::String  = "",

        tipjointtag::String  = "",
        tipjoint   ::Symbol = :none,
        id  ::Int = -1,
    )
        if path.closed || embedded
            tipjoint = :none
        end
        return new(id, path, closed, embedded, shape, tag, jointtag, tipjointtag, tipjoint)
    end
end