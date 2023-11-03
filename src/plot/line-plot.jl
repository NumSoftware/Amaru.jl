# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

const _line_style_list = [:solid, :dash, :dashdot]
const _marker_list = [:none, :circle, :square, :triangle, :utriangle, :cross, :xcross, :diamond, :pentagon, :hexagon, :star]


mutable struct LinePlot<:DataSeriesPlot
    X     ::Array
    Y     ::Array
    ls    ::Symbol
    lw    ::Float64
    lc    ::Union{Symbol,Tuple,Nothing,Vec3}
    marker::Symbol
    ms    ::Float64
    label ::String
    dash  ::Vector{Float64}
    # xmult::String
    # ymult::String

    function LinePlot( X::AbstractArray, Y::AbstractArray; args...)

        args = checkargs( args, 
            ArgInfo( (:ls, :linestyle), "Line style", default=:solid, values=_line_style_list ),
            ArgInfo( (:lc, :linecolor,:color), "Line lc", default=:default, values=[:default; _colors_list]),
            ArgInfo( (:lw, :lineweight), "Line weight", default=0.5,  condition=:(lw>0) ),
            ArgInfo( :marker, "Marker shape", default=:none,  values=_marker_list ),
            ArgInfo( (:ms, :markersize), "Marker size", default=3, condition=:(ms>0) ),
            ArgInfo( (:mfc, :markerfacecolor), "Marker face lc", default=:default ),
            ArgInfo( :label, "Legend label", default=""),
            # ArgInfo( :xlims, "Limits for the x axis", default=:default ),
            # ArgInfo( :ylims, "Limits for the y axis", default=:default ),
        )

        this = new(X, Y, args.ls, args.lw, args.lc, args.marker, args.ms, args.label)
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

    if p.lc===:default
        p.lc = _colors_dict[_default_colors[chart.icolor]]
        chart.icolor += 1
    else
        p.lc = _colors_dict[p.lc]
    end

    set_matrix(cc, CairoMatrix([1, 0, 0, 1, 0, 0]...)) 
    set_source_rgb(cc, p.lc...)
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

        draw_marker(cc, x, y, p.ms, p.marker)
    end
end


function draw_polygon(cc::CairoContext, x, y, n, l; angle=0)
    Δθ = 360/n
    minθ = angle + 90
    maxθ = angle + 360 + 90

    for θ in minθ:Δθ:maxθ
        xi = x + l*cosd(θ)
        yi = y - l*sind(θ)
        if θ==angle
            move_to(cc, xi, yi)
        else
            line_to(cc, xi, yi)
        end
    end

    close_path(cc)
    fill(cc)
end


function draw_star(cc::CairoContext, x, y, n, l; angle=0)
    Δθ = 360/n/2
    minθ = angle + 90
    maxθ = angle + 360 + 90


    for (i,θ) in enumerate(minθ:Δθ:maxθ)
        if i%2==1
            xi = x + l*cosd(θ)
            yi = y - l*sind(θ)
        else
            xi = x + 0.5*l*cosd(θ)
            yi = y - 0.5*l*sind(θ)
        end
        if θ==angle
            move_to(cc, xi, yi)
        else
            line_to(cc, xi, yi)
        end
    end

    close_path(cc)
    fill(cc)
end


function draw_marker(cc::CairoContext, x, y, size, marker)
    radius = size/2

    if marker==:circle
        arc(cc, x, y, radius, 0, 2*pi)
        fill(cc)
    elseif marker==:square
        draw_polygon(cc, x, y, 4, 1.2*radius, angle=45)
    elseif marker==:diamond
        draw_polygon(cc, x, y, 4, 1.2*radius, angle=0)
    elseif marker==:triangle
        draw_polygon(cc, x, y, 3, 1.3*radius, angle=0)
    elseif marker==:utriangle
        draw_polygon(cc, x, y, 3, 1.3*radius, angle=180)
    elseif marker==:pentagon
        draw_polygon(cc, x, y, 5, 1.1*radius, angle=0)
    elseif marker==:hexagon
        draw_polygon(cc, x, y, 5, 1.1*radius, angle=0)
    elseif marker==:star
        draw_star(cc, x, y, 5, 1.25*radius, angle=0)
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