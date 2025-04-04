export PathFunction

mutable struct PathFunCmd
    key::Symbol
    points::Array{Vec2}
end

function (cmd::PathFunCmd)(t::Float64) 
    if cmd.key==:M
        return cmd.points[1]
    elseif cmd.key==:L
        return cmd.points[1] + t*(cmd.points[2]-cmd.points[1])
    elseif cmd.key==:Q
        return (1-t)^2*cmd.points[1] + 2*(1-t)*t*cmd.points[2] + t^2*cmd.points[3]
    elseif cmd.key==:C
        return (1-t)^3*cmd.points[1] + 3*(1-t)^2*t*cmd.points[2] + 3*(1-t)*t^2*cmd.points[3] + t^3*cmd.points[4]
    else
        error("PathFunCmd: Invalid command key $(cmd.key)")
    end
end


# function to find the parameter t for a given x
function find_t(cmd::PathFunCmd, x::Real)
    if cmd.key==:M
        return 0.0
    elseif cmd.key==:L
        x1 = cmd.points[1][1]
        x2 = cmd.points[2][1]

        return (x-x1)/(x2-x1)
    elseif cmd.key==:Q
        x1 = cmd.points[1][1]
        x2 = cmd.points[2][1]
        x3 = cmd.points[3][1]

        a = x1 - 2*x2 + x3
        b = -2*x1 + 2*x2
        c = x1 - x

        if a==0
            t = -c/b
        else
            disc = b^2 - 4*a*c
            t = (-b + √disc)/(2*a)

            if t<0 || t>1
                t = (-b - √disc)/(2*a)
            end
        end

        return t
    elseif cmd.key==:C
        x1 = cmd.points[1][1]
        x2 = cmd.points[2][1]
        x3 = cmd.points[3][1]
        x4 = cmd.points[4][1]


        f(t) = (1-t)^3*x1 + 3*(1-t)^2*t*x2 + 3*(1-t)*t^2*x3 + t^3*x4 - x
        df(t) = 3*(1-t)^2*(x2-x1) + 6*(1-t)*t*(x3-x2) + 3*t^2*(x4-x3)
    
        tol = 0.001
        maxits = 10
        t = (x-x1)/(x4-x1)
        for i in 1:maxits
            t = t - f(t)/df(t)
            abs(f(t))<tol && break
            i==10 && error("Failed to converge")
        end

        return t
    else
        error("PathFunCmd: Invalid command key $(cmd.key)")
    end
end

function interpolate(cmd::PathFunCmd, x::Real)
    if cmd.key==:M
        return cmd.points[1][2]
    elseif cmd.key==:L
        x1, y1 = cmd.points[1]
        x2, y2 = cmd.points[2]
        x<=x1 && return y1
        x>=x2 && return y2

        t = (x-x1)/(x2-x1)
        return y1 + t*(y2-y1)
    elseif cmd.key==:Q
        x1, y1 = cmd.points[1]
        x2, y2 = cmd.points[2]
        x3, y3 = cmd.points[3]
        x<=x1 && return y1
        x>=x3 && return y3

        t = find_t(cmd, x)

        return (1-t)^2*y1 + 2*(1-t)*t*y2 + t^2*y3
        
    elseif cmd.key==:C
        x1, y1 = cmd.points[1]
        x2, y2 = cmd.points[2]
        x3, y3 = cmd.points[3]
        x4, y4 = cmd.points[4]
        x<=x1 && return y1
        x>=x4 && return y4

        t = find_t(cmd, x)

        return (1-t)^3*y1 + 3*(1-t)^2*t*y2 + 3*(1-t)*t^2*y3 + t^3*y4
    else
        error("PathFunCmd: Invalid command key $(cmd.key)")
    end
end

function derive(cmd::PathFunCmd, x::Real)
    xini = cmd.points[1][1]
    xend = cmd.points[end][1]
    x<xini && return 0.0
    x>xend && return 0.0

    if cmd.key==:M
        return 0.0
    elseif cmd.key==:L
        x1, y1 = cmd.points[1]
        x2, y2 = cmd.points[2]

        return (y2-y1)/(x2-x1)
    elseif cmd.key==:Q
        x1, y1 = cmd.points[1]
        x2, y2 = cmd.points[2]
        x3, y3 = cmd.points[3]

        t = find_t(cmd, x)
        dydt = 2*(1-t)*(y2-y1) + 2*t*(y3-y2)
        dxdt = 2*(1-t)*(x2-x1) + 2*t*(x3-x2)
        
        return dydt/dxdt
    elseif cmd.key==:C
        x1, y1 = cmd.points[1]
        x2, y2 = cmd.points[2]
        x3, y3 = cmd.points[3]
        x4, y4 = cmd.points[4]

        t = find_t(cmd, x)
        dydt = 3*(1-t)^2*(y2-y1) + 6*(1-t)*t*(y3-y2) + 3*t^2*(y4-y3)
        dxdt = 3*(1-t)^2*(x2-x1) + 6*(1-t)*t*(x3-x2) + 3*t^2*(x4-x3)
        
        return dydt/dxdt
    else
        error("Invalid command key $(cmd.key)")
    end
end


## Path function ##

struct PathFunction
    points::Vector{Vec2}
    cmds::Vector{PathFunCmd}

    function PathFunction()
        return new(Vector{Vec2}(), Vector{PathFunCmd}())
    end

    function PathFunction(list::Union{Symbol, Real}...)
        this = new(Vector{Vec2}(), Vector{PathFunCmd}())
        append!(this, list...)
        return this
    end

end


function Base.append!(f::PathFunction, list::Union{Symbol, Real}...)
    n   = length(list)
    idx = 1

    while idx<=n 
        key = list[idx]
        if key==:M
            n>=idx+2 || error("PathFunction: M command requires at least two numbers")
            p1 = Vec2(list[idx+1], list[idx+2])
            push!(f.points, p1)
            push!(f.cmds, PathFunCmd(key, [p1]))
            idx += 3
        elseif key==:L
            n>=idx+2 || error("PathFunction: L command requires at least two numbers")
            p1 = f.points[end]
            p2 = Vec2(list[idx+1], list[idx+2])
            push!(f.points, p2)
            push!(f.cmds, PathFunCmd(key, [p1, p2]))
            idx += 3
        elseif key==:Q
            n>=idx+4 || error("PathFunction: Q command requires at least four numbers")
            p1 = f.points[end]
            p2 = Vec2(list[idx+1], list[idx+2])
            p3 = Vec2(list[idx+3], list[idx+4])
            push!(f.points, p3)
            push!(f.cmds, PathFunCmd(key, [p1, p2, p3]))
            idx += 5
        elseif key==:C
            n>=idx+6 || error("PathFunction: C command requires at least six numbers")
            p1 = f.points[end]
            p2 = Vec2(list[idx+1], list[idx+2])
            p3 = Vec2(list[idx+3], list[idx+4])
            p4 = Vec2(list[idx+5], list[idx+6])
            push!(f.points, p4)
            push!(f.cmds, PathFunCmd(key, [p1, p2, p3, p4]))
            idx += 7
        else
            error("PathFunction: Invalid command $key")
        end
    end

    return f
end


function Base.extrema(pathf::PathFunction)
    x0 = pathf.points[1][1]
    xend = pathf.points[end][1]
    X = range(x0, xend, length=100)
    Y = pathf.(X)
    return extrema(Y)
end


function interpolate(path::PathFunction, x::Real)
    @assert !isnan(x)
    if x<=path.points[1][1]
        return path.points[1][2]
    end

    if x>=path.points[end][1]
        return path.points[end][2]
    end

    # search for the cmd that contains x
    cmd = path.cmds[1]
    for c in path.cmds[2:end]
        if x<=c.points[end][1]
            cmd = c
            break
        end
    end

    # interpolate 
    return interpolate(cmd, x)

end

@inline function (path::PathFunction)(x::Real)
    return interpolate(path, x)
end


function derive(path::PathFunction, x::Real)

    if x<path.points[1][1] || x>path.points[end][1]
        return 0.0
    end

    # search for the cmd that contains x
    cmd = path.cmds[1]
    for c in path.cmds[2:end]
        if x<=c.points[end][1]
        # if x<=lastpoint(c)[1]
            cmd = c
            break
        end
    end

    # interpolate 
    return derive(cmd, x)

end
    

# F = PathFunction(:M, 0, 0, :Q, 0.4, 1, 1, 1 )
# F = PathFunction(:M, 0, 0, :Q, 0.5, 1, 1, 1, :C, 1.3, 1, 1.5, 0, 2, 0, :L, 3, 1 )

# @show interpolate(F, 2.5)
# @show derive(F, 2.5)


# using Amaru

# chart = Chart()

# X = 0:0.01:4
# Y = [ interpolate(F, x) for x in X ]
# D = [ derive(F, x) for x in X ]

# series = [ 
#     LineSeries(X, Y)
#     LineSeries(X, D)
#  ]
# addseries!(chart, series)
# save(chart, "chart.pdf")