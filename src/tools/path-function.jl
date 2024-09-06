

export PathFunction

abstract type PathFunCmd
end


# Move command
# ============

struct MovePFCmd<:PathFunCmd
    p1::Vec2
end

lastpoint(cmd::MovePFCmd) = cmd.p1
(mc::MovePFCmd)(t::Float64) = mc.p1.coord


# Line command
# ============

struct LinePFCmd<:PathFunCmd
    p1::Vec2
    p2::Vec2
end

lastpoint(cmd::LinePFCmd) = cmd.p2
(c::LinePFCmd)(t::Float64) = c.p1 + t*(c.p2-c.p1) 

function interpolate(c::LinePFCmd, x::Real)
    x1, y1 = c.p1
    x2, y2 = c.p2
    x<=x1 && return y1
    x>=x2 && return y2

    t = (x-x1)/(x2-x1)
    x2==x1 && return y1
    return y1 + t*(y2-y1)
end

function derive(c::LinePFCmd, x::Real)
    x1, y1 = c.p1
    x2, y2 = c.p2
    x<x1 && return 0.0
    x>x2 && return 0.0

    return (y2-y1)/(x2-x1)
end


# Quadratic command
# =================

struct QuadraticPFCmd<:PathFunCmd
    p1::Vec2
    p2::Vec2
    p3::Vec2
end

lastpoint(cmd::QuadraticPFCmd) = cmd.p3
(c::QuadraticPFCmd)(t::Float64) = (1-t)^2*c.p1 + 2*(1-t)*t*c.p2 + t^2*c.p3

function find_t(c::QuadraticPFCmd, x::Real)
    x1 = c.p1[1]
    x2 = c.p2[1]
    x3 = c.p3[1]

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
end

function interpolate(c::QuadraticPFCmd, x::Real)
    x1, y1 = c.p1
    x2, y2 = c.p2
    x3, y3 = c.p3
    x<=x1 && return y1
    x>=x3 && return y3

    t = find_t(c, x)

    return (1-t)^2*y1 + 2*(1-t)*t*y2 + t^2*y3
end

function derive(c::QuadraticPFCmd, x::Real)
    x1, y1 = c.p1
    x2, y2 = c.p2
    x3, y3 = c.p3
    x<x1 && return 0.0
    x>x3 && return 0.0

    t = find_t(c, x)
    dydt = 2*(1-t)*(y2-y1) + 2*t*(y3-y2)
    dxdt = 2*(1-t)*(x2-x1) + 2*t*(x3-x2)
    
    return dydt/dxdt
end


# Cubic Bezier command
# ====================

struct BezierPFCmd<:PathFunCmd
    p1::Vec2
    p2::Vec2
    p3::Vec2
    p4::Vec2
end

lastpoint(cmd::BezierPFCmd) = cmd.p4
(c::BezierPFCmd)(t::Float64) = (1-t)^3*c.p1 + 3*(1-t)^2*t*c.p2 + 3*(1-t)*t^2*c.p3 + t^3*c.p4


function find_t(c::BezierPFCmd, x::Real)
    x1 = c.p1[1]
    x2 = c.p2[1]
    x3 = c.p3[1]
    x4 = c.p4[1]

    f(t) = (1-t)^3*x1 + 3*(1-t)^2*t*x2 + 3*(1-t)*t^2*x3 + t^3*x4 - x
    df(t) = 3*(1-t)^2*(x2-x1) + 6*(1-t)*t*(x3-x2) + 3*t^2*(x4-x3)

    tol = 0.001
    maxits = 10
    t = (x-x1)/(x4-x1)
    # t = 0.5
    for i in 1:maxits
        t = t - f(t)/df(t)
        abs(f(t))<tol && break
        i==10 && error("Failed to converge")
    end

    return t
end

function interpolate(c::BezierPFCmd, x::Real)
    x1, y1 = c.p1
    x2, y2 = c.p2
    x3, y3 = c.p3
    x4, y4 = c.p4
    x<=x1 && return y1
    x>=x4 && return y4

    t = find_t(c, x)

    return (1-t)^3*y1 + 3*(1-t)^2*t*y2 + 3*(1-t)*t^2*y3 + t^3*y4
end

function derive(c::BezierPFCmd, x::Real)
    x1, y1 = c.p1
    x2, y2 = c.p2
    x3, y3 = c.p3
    x4, y4 = c.p4
    x<=x1 && return 0.0
    x>=x4 && return 0.0

    t = find_t(c, x)
    dydt = 3*(1-t)^2*(y2-y1) + 6*(1-t)*t*(y3-y2) + 3*t^2*(y4-y3)
    dxdt = 3*(1-t)^2*(x2-x1) + 6*(1-t)*t*(x3-x2) + 3*t^2*(x4-x3)
    
    return dydt/dxdt
end


# Path function
# =============

struct PathFunction
    points::Vector{Vec2}
    cmds::Vector{PathFunCmd}

    function PathFunction()
        return new(Vector{Vec2}(), Vector{PathFunCmd}())
    end
end

lastpoint(path::PathFunction) = lastpoint(path.cmds[end])


# function to create a path from a list of commands
function PathFunction(list::Union{Symbol, Real}...)
    f = PathFunction()
    append!(f, list...)
    return f
end


function Base.append!(f::PathFunction, list::Union{Symbol, Real}...)
    n   = length(list)
    idx = 1

    while idx<=n 
        item = list[idx]
        if item==:M
            n>=idx+2 || error("PathFunction: M command requires at least two numbers")
            p1 = Vec2(list[idx+1], list[idx+2])
            push!(f.points, p1)
            push!(f.cmds, MovePFCmd(p1))
            idx += 3
        elseif item==:L
            n>=idx+2 || error("PathFunction: L command requires at least two numbers")
            p1 = f.points[end]
            p2 = Vec2(list[idx+1], list[idx+2])
            push!(f.points, p2)
            push!(f.cmds, LinePFCmd(p1, p2))
            idx += 3
        elseif item==:Q
            n>=idx+4 || error("PathFunction: Q command requires at least four numbers")
            p1 = f.points[end]
            p2 = Vec2(list[idx+1], list[idx+2])
            p3 = Vec2(list[idx+3], list[idx+4])
            push!(f.points, p3)
            push!(f.cmds, QuadraticPFCmd(p1, p2, p3))
            idx += 5
        elseif item==:C
            n>=idx+6 || error("PathFunction: C command requires at least six numbers")
            p1 = f.points[end]
            p2 = Vec2(list[idx+1], list[idx+2])
            p3 = Vec2(list[idx+3], list[idx+4])
            p4 = Vec2(list[idx+5], list[idx+6])
            push!(f.points, p4)
            push!(f.cmds, BezierPFCmd(p1, p2, p3, p4))
            idx += 7
        else
            error("PathFunction: Invalid command $item")
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
    if x<=path.cmds[1].p1[1]
        return path.cmds[1].p1[2]
    end

    lastp = lastpoint(path.cmds[end])
    if x>=lastp[1]
        return lastp[2]
    end

    # search for the cmd that contains x
    cmd = path.cmds[1]
    for c in path.cmds[2:end]
        if x<=lastpoint(c)[1]
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
    if x<path.cmds[1].p1[1]
        return 0.0
    end

    lastp = lastpoint(path.cmds[end])
    if x>lastp[1]
        return 0.0
    end

    # search for the cmd that contains x
    cmd = path.cmds[1]
    for c in path.cmds[2:end]
        if x<=lastpoint(c)[1]
            cmd = c
            break
        end
    end

    # interpolate 
    return derive(cmd, x)

end
    

# F = PathFunction([:M, 0, 0, :Q, 0.4, 1, 1, 1] )
# F = PathFunction([:M, 0, 0, :Q, 0.5, 1, 1, 1, :B, 1.3, 1, 1.5, 0, 2, 0, :L, 3, 1] )

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