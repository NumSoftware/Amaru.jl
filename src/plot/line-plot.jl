# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

const _line_style_list = [:none, :solid, :dot, :dash, :dashdot]
const _marker_list = [:none, :circle, :square, :triangle, :utriangle, :cross, :xcross, :diamond, :pentagon, :hexagon, :star]


mutable struct LinePlot<:DataSeriesPlot
    X     ::Array
    Y     ::Array
    ls    ::Symbol
    lw    ::Float64
    linecolor::Union{Symbol,Tuple}
    marker::Symbol
    markersize::Float64
    markercolor::Union{Symbol,Tuple}
    markerstrokecolor::Union{Symbol,Tuple}
    label ::String
    dash  ::Vector{Float64}
    order::Union{Int,Nothing}

    function LinePlot( X::AbstractArray, Y::AbstractArray; args...)

        args = checkargs(args, func_params(LinePlot), aliens=false)

        if args.linecolor!==:default
            linecolor = get_color(args.linecolor)
            markercolor = get_color(args.markercolor, linecolor)
            markerstrokecolor = get_color(args.markerstrokecolor, linecolor)
        else
            linecolor = args.linecolor
            markercolor = args.markercolor
            markerstrokecolor = args.markerstrokecolor
        end
        n = min(length(X), length(Y))

        ls = length(args.dash)>0 ? :dash : args.ls

        this = new(X[1:n], Y[1:n], ls, args.lw, linecolor, args.marker, args.markersize, markercolor, markerstrokecolor, args.label, args.dash, args.order)
        return this
    end
end

func_params(::Type{LinePlot}) = [
    FunInfo( :LinePlot, "Creates a customizable `LinePlot` instance.", (Array, Array)),
    ArgInfo( (:ls, :linestyle), "Line style", :solid, values=_line_style_list ),
    ArgInfo( :dash, "Dash pattern", Float64[] ),
    ArgInfo( (:linecolor, :lc, :color), "Line linecolor", :default),
    ArgInfo( (:lw, :lineweight), "Line weight", 0.5,  condition=:(lw>0) ),
    ArgInfo( :marker, "Marker shape", :none,  values=_marker_list ),
    ArgInfo( (:markersize, :ms), "Marker size", 2.5, condition=:(markersize>0) ),
    ArgInfo( (:markercolor, :mc), "Marker color", :white ),
    ArgInfo( (:markerstrokecolor, :msc), "Marker stroke color", :default ),
    ArgInfo( :label, "Legend label", ""),
    ArgInfo( :order, "Order fo drawing", nothing),
]
@doc make_doc(LinePlot) LinePlot()

function data2user(c::Chart, x, y)
    Xmin, Ymin, Xmax, Ymax = c.canvas.box
    xmin, ymin, xmax, ymax = c.canvas.limits
    xD = Xmin + (Xmax-Xmin)/(xmax-xmin)*(x-xmin)
    yD = Ymin + (Ymax-Ymin)/(ymax-ymin)*(ymax-y)
    return xD, yD
end

function configure!(chart::Chart, p::LinePlot)
    if length(p.dash)==0
        if p.ls==:dash
            p.dash = [4.0, 2.4]*p.lw
        elseif p.ls==:dashdot
            p.dash = [2.0, 1.0, 2.0, 1.0]*p.lw
        elseif p.ls==:dot
            p.dash = [1.0, 1.0]*p.lw
        end
    end
end

function draw!(chart::Chart, cc::CairoContext, p::LinePlot)

    p.markercolor = get_color(p.markercolor, p.linecolor)
    p.markerstrokecolor = get_color(p.markerstrokecolor, p.linecolor)

    set_matrix(cc, CairoMatrix([1, 0, 0, 1, 0, 0]...)) 
    set_source_rgb(cc, p.linecolor...)
    set_line_width(cc, p.lw)

    # Plot lines
    n = length(p.X)
    X = p.X*chart.xaxis.mult
    Y = p.Y*chart.yaxis.mult
    x1, y1 = data2user(chart, X[1], Y[1])


    if p.ls!==:none
        p.ls!=:solid && set_dash(cc, p.dash)

        for i in 2:n
            x, y = data2user(chart, X[i], Y[i])
            move_to(cc, x1, y1); line_to(cc, x, y); stroke(cc)
            x1, y1 = x, y
        end
    end

    set_dash(cc, Float64[])
    
    # Plot markers
    for (x,y) in zip(X, Y)
        x, y = data2user(chart, x, y)
        draw_marker(cc, x, y, p.marker, p.markersize, p.markercolor, p.markerstrokecolor)
    end
end


function draw_polygon(cc::CairoContext, x, y, n, length, color, strokecolor; angle=0)
    Δθ = 360/n
    minθ = angle + 90
    maxθ = angle + 360 + 90

    for θ in minθ:Δθ:maxθ
        xi = x + length*cosd(θ)
        yi = y - length*sind(θ)
        if θ==angle
            move_to(cc, xi, yi)
        else
            line_to(cc, xi, yi)
        end
    end

    close_path(cc)
    set_source_rgb(cc, color...)
    fill_preserve(cc)
    set_source_rgb(cc, strokecolor...)
    stroke(cc)
end


function draw_star(cc::CairoContext, x, y, n, length, color, strokecolor; angle=0)
    Δθ = 360/n/2
    minθ = angle + 90
    maxθ = angle + 360 + 90


    for (i,θ) in enumerate(minθ:Δθ:maxθ)
        if i%2==1
            xi = x + length*cosd(θ)
            yi = y - length*sind(θ)
        else
            xi = x + 0.5*length*cosd(θ)
            yi = y - 0.5*length*sind(θ)
        end
        if θ==angle
            move_to(cc, xi, yi)
        else
            line_to(cc, xi, yi)
        end
    end

    close_path(cc)
    set_source_rgb(cc, color...)
    fill_preserve(cc)
    set_source_rgb(cc, strokecolor...)
    stroke(cc)
end


function draw_marker(cc::CairoContext, x, y, marker, size, color, strokecolor)
    radius = size/2

    if marker==:circle
        arc(cc, x, y, radius, 0, 2*pi)
        set_source_rgb(cc, color...)
        fill_preserve(cc)
        set_source_rgb(cc, strokecolor...)
        stroke(cc)
    elseif marker==:square
        draw_polygon(cc, x, y, 4, 1.2*radius, color, strokecolor, angle=45)
    elseif marker==:diamond
        draw_polygon(cc, x, y, 4, 1.2*radius, color, strokecolor, angle=0)
    elseif marker==:triangle
        draw_polygon(cc, x, y, 3, 1.3*radius, color, strokecolor, angle=0)
    elseif marker==:utriangle
        draw_polygon(cc, x, y, 3, 1.3*radius, color, strokecolor, angle=180)
    elseif marker==:pentagon
        draw_polygon(cc, x, y, 5, 1.1*radius, color, strokecolor, angle=0)
    elseif marker==:hexagon
        draw_polygon(cc, x, y, 6, 1.1*radius, color, strokecolor, angle=0)
    elseif marker==:star
        draw_star(cc, x, y, 5, 1.25*radius, color, strokecolor, angle=0)
    elseif marker==:cross
        radius = 1.35*radius
        set_line_width(cc, radius/3)
        move_to(cc, x, y-radius)
        line_to(cc, x, y+radius)
        stroke(cc)
        move_to(cc, x-radius, y)
        line_to(cc, x+radius, y)
        stroke(cc)
    elseif marker==:xcross
        radius = 1.35*radius
        set_line_width(cc, radius/3)
        move_to(cc, x-radius, y-radius)
        line_to(cc, x+radius, y+radius)
        stroke(cc)
        move_to(cc, x+radius, y-radius)
        line_to(cc, x-radius, y+radius)
        stroke(cc)
    end
end