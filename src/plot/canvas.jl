# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct Canvas<:FigureComponent
    width::Float64
    height::Float64
    limits::Vector{Float64}  # [ xmin, ymin, xmax, ymax ] in data coordinates
    box::Vector{Float64}     # [ x0, y0, x1, y1 ] in svg coordinates
    function Canvas()
        return new(0.0, 0.0)
    end
end

function data2user(canvas::Canvas, x::Float64, y::Float64)
    Xmin, Ymin, Xmax, Ymax = canvas.box
    xmin, ymin, xmax, ymax = canvas.limits
    X = Xmin + (Xmax-Xmin)/(xmax-xmin)*(x-xmin)
    Y = Ymin + (Ymax-Ymin)/(ymax-ymin)*(ymax-y)
    return X, Y
end

function user2data(canvas::Canvas, X::Float64, Y::Float64)
    Xmin, Ymin, Xmax, Ymax = canvas.box
    xmin, ymin, xmax, ymax = canvas.limits
    x = xmin + (xmax-xmin)/(Xmax-Xmin)*(X-Xmin)
    y = ymax + (ymax-ymin)/(Ymax-Ymin)*(Ymin-Y)
    return x, y
end