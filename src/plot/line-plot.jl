# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

const _line_style_list = [:solid, :dash, :dashdot]
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
    # xmult::String
    # ymult::String

    function LinePlot( X::AbstractArray, Y::AbstractArray; args...)

        args = checkargs( args, 
            [
                ArgInfo( (:ls, :linestyle), "Line style", default=:solid, values=_line_style_list ),
                ArgInfo( (:linecolor, :lc, :color), "Line linecolor", default=:default),
                ArgInfo( (:lw, :lineweight), "Line weight", default=0.5,  condition=:(lw>0) ),
                ArgInfo( :marker, "Marker shape", default=:none,  values=_marker_list ),
                ArgInfo( (:markersize, :ms), "Marker size", default=2.5, condition=:(markersize>0) ),
                ArgInfo( (:markercolor, :mc), "Marker color", default=:white ),
                ArgInfo( (:markerstrokecolor, :msc), "Marker stroke color", default=:default ),
                ArgInfo( :label, "Legend label", default=""),
                # ArgInfo( :xlims, "Limits for the x axis", default=:default ),
                # ArgInfo( :ylims, "Limits for the y axis", default=:default ),
            ],
            checkwrong=true
        )

        if args.linecolor!==:default
            linecolor = get_color(args.linecolor)
            markercolor = get_color(args.markercolor, linecolor)
            markerstrokecolor = get_color(args.markerstrokecolor, linecolor)
        else
            linecolor = args.linecolor
            markercolor = args.markercolor
            markerstrokecolor = args.markerstrokecolor
        end

        this = new(X, Y, args.ls, args.lw, linecolor, args.marker, args.markersize, markercolor, markerstrokecolor, args.label)
        return this
    end
end


function data2user(c::Chart, x, y)
    Xmin, Ymin, Xmax, Ymax = c.canvas.box
    xmin, ymin, xmax, ymax = c.canvas.limits
    xD = Xmin + (Xmax-Xmin)/(xmax-xmin)*(x-xmin)
    yD = Ymin + (Ymax-Ymin)/(ymax-ymin)*(ymax-y)
    return xD, yD
end

function configure!(chart::Chart, p::LinePlot)
    if p.ls==:dash
        reflen = 0.01*minimum(chart.figsize)
        p.dash = [reflen, 0.7*reflen]
    else
        p.dash = []
    end
end

function draw!(chart::Chart, cc::CairoContext, p::LinePlot)

    if p.linecolor===:default # update colors
        p.linecolor = _colors_dict[_default_colors[chart.icolor]]
        p.markercolor = get_color(p.markercolor, p.linecolor)
        p.markerstrokecolor = get_color(p.markerstrokecolor, p.linecolor)
        chart.icolor += 1
    end

    set_matrix(cc, CairoMatrix([1, 0, 0, 1, 0, 0]...)) 
    set_source_rgb(cc, p.linecolor...)
    set_line_width(cc, p.lw)
    if p.ls!=:solid
        set_dash(cc, p.dash)
    end

    # Plot lines
    n = length(p.X)
    X = p.X*chart.xaxis.mult
    Y = p.Y*chart.yaxis.mult
    x1, y1 = data2user(chart, X[1], Y[1])

    for i in 2:n
        x, y = data2user(chart, X[i], Y[i])
        move_to(cc, x1, y1); line_to(cc, x, y); stroke(cc)
        x1, y1 = x, y
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
    end
end