export Path, addinset!

mutable struct PathCmd
    key::Symbol
    points::Array{Point}
end

function Base.copy(cmd::PathCmd)
    return PathCmd(cmd.key, copy.(cmd.points))
end


function Base.length(cmd::PathCmd)
    if cmd.key==:M
        return 0.0
    elseif cmd.key==:L
        return norm(cmd.points[2].coord-cmd.points[1].coord)
    elseif cmd.key==:A
        X1 = cmd.points[1].coord
        X2 = cmd.points[2].coord
        X3 = cmd.points[3].coord
        r = norm(X2-X1)
        θ = acos(clamp(dot(X1-X2, X3-X2)/r^2, -1, 1))
        len = r*θ
        return len
    elseif cmd.key==:C
        X1 = cmd.points[1].coord
        X2 = cmd.points[2].coord
        X3 = cmd.points[3].coord
        X4 = cmd.points[4].coord

        deriv(t) = 3*(1-t)^2*(X2-X1) + 6*(1-t)*t*(X3-X2) + 3*t^2*(X4-X3)
        ipoints = [0.5 - √(3/5)/2, 0.5, 0.5 + √(3/5)/2]
        weights = [5/18, 8/18, 5/18]
        return sum([weights[i]*norm(deriv(ipoints[i])) for i in 1:3])
    else
        error("Invalid command key $(cmd.key)")
    end

end


function (cmd::PathCmd)(t::Float64)
    if cmd.key==:M
        return cmd.points[1].coord
    elseif cmd.key==:L
        return cmd.points[1].coord + t*(cmd.points[2].coord-cmd.points[1].coord)
    elseif cmd.key==:A
        X1 = cmd.points[1].coord
        X2 = cmd.points[2].coord
        X3 = cmd.points[3].coord
        axis = cross(X1-X2, X3-X2)
        axis = axis/norm(axis)
        r = norm(X2-X1)
        θ = t*acos(clamp(dot(X1-X2, X3-X2)/r^2, -1, 1))
        R = Quaternion(cos(θ/2), axis[1]*sin(θ/2), axis[2]*sin(θ/2), axis[3]*sin(θ/2))
        P = X1 + R*(X2-X1)*conj(R)
        return P
    elseif cmd.key==:C
        X1 = cmd.points[1].coord
        X2 = cmd.points[2].coord
        X3 = cmd.points[3].coord
        X4 = cmd.points[4].coord
        (1-t)^3*X1 + 3*(1-t)^2*t*X2 + 3*(1-t)*t^2*X3 + t^3*X4
    else
        error("Invalid command key $(cmd.key)")
    end
end


mutable struct Path
    points::Array{Point}
    cmds::Array{PathCmd}
    closed::Bool
    len::Array{Float64,1} # normalized cumulative length

    function Path(points::Array{Point}, cmds::Array{PathCmd}, closed::Bool, len::Array{Float64,1})
        return new(points, cmds, closed, len)
    end
end


function Base.copy(path::Path)
    # points = copy.(path.points)
    cmds = copy.(path.cmds) # also copies the points
    points = Point[]
    for cmd in cmds
        if cmd.key==:M
            push!(points, cmd.points[1])
        else 
            cmd.points[1] = points[end]
            append!(points, cmd.points[2:end])
        end
    end

    return Path(points, cmds, path.closed, copy(path.len))

end


# todo: check this function.. currently not used
function path_from_numbers(tokens...; closed=false)
    n  = length(tokens)
    idx = 1
    cmds = PathCmd[]
    local startpoint, endpoint

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
    # any( t->isa(t, Point), tokens ) || return path_from_numbers(tokens...; closed=closed)

    n  = length(tokens)
    idx = 1
    cmds = PathCmd[]
    points = Point[]

    # local startpoint, endpoint

    while idx<=n
        token = tokens[idx]
        if token==:M
            p1 = tokens[idx+1]
            cmd = PathCmd(token, [ p1 ])
            # cmd = MoveCmd(p1)
            push!(points, p1)
            push!(cmds, cmd)
            idx+=2
        elseif token==:L
            p1 = points[end]
            p2 = tokens[idx+1]
            cmd = PathCmd(token, [ p1, p2 ])
            # cmd = LineCmd(p1, p2)
            push!(points, p2)
            push!(cmds, cmd)
            idx+=2
        elseif token==:A
            p1 = points[end]
            p2 = tokens[idx+1]
            p3 = tokens[idx+2]
            cmd = PathCmd(token, [p1, p2, p3])
            # cmd = ArcCmd(p1, p2, p3)
            push!(points, p2)
            push!(points, p3)
            push!(cmds, cmd)
            idx+=3
        elseif token==:C
            p1 = points[end]
            p2 = tokens[idx+1]
            p3 = tokens[idx+2]
            p4 = tokens[idx+3]
            cmd = PathCmd(token, [p1, p2, p3, p4])
            # cmd = BezierCmd(p1, p2, p3, p4)
            push!(points, p2)
            push!(points, p3)
            push!(points, p4)
            push!(cmds, cmd)
            idx+=4
        else
            error("Invalid token $token")
        end
        
    end

    if closed && points[1]!==points[end]
        cmd = LineCmd(endpoint, startpoint)
        push!(cmds, cmd)
    end

    if !closed && points[1]===points[end]
        closed = true
    end

    # normalized length
    L = cumsum( length(pc) for pc in cmds )
    L = L./L[end] # normalize

    return Path(points, cmds, closed, L)
end


function (path::Path)(t::Float64) 
    i = searchsortedfirst(path.len, t)
    if i==1
        return path.cmds[1](0.0)
    elseif i>length(path.len)
        return path.cmds[end](1.0)
    else
        tt = (t-path.len[i-1])/(path.len[i]-path.len[i-1])
        return path.cmds[i](tt)
    end
end


mutable struct SubPath<:GeoEntity
    id      ::Int
    path    ::Path
    closed  ::Bool
    embedded::Bool
    shape   ::CellShape
    tag     ::String
    jointtag::String
    tiptag  ::String
    tips    ::Symbol # :none, :start, :end, :both

    function SubPath(
        path::Path;
        closed::Bool  = false,
        embedded::Bool  = false,
        shape   ::CellShape = LIN3,
        tag     ::String  = "",
        jointtag::String  = "",

        tiptag::String  = "",
        tips::Symbol = :none,
        id  ::Int = -1,
    )
        if path.closed || embedded
            tips = :none
        end
        return new(id, path, closed, embedded, shape, tag, jointtag, tiptag, tips)
    end
end


function Base.copy(sp::SubPath)
    return SubPath(copy(sp.path), closed=sp.closed, embedded=sp.embedded, shape=sp.shape, tag=sp.tag, jointtag=sp.jointtag, tiptag=sp.tiptag, tips=sp.tips)
end