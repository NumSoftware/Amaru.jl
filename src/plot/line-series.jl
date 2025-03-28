# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

const _line_style_list = [:none, :solid, :dot, :dash, :dashdot]
const _marker_list = [:none, :circle, :square, :triangle, :utriangle, :cross, :xcross, :diamond, :pentagon, :hexagon, :star]


LineSeries_params = [
    FunInfo(:LineSeries, "Creates a customizable `LineSeries` instance to be used in a `Chart`"),
    ArgInfo(:X, "Array of x-coordinates"),
    ArgInfo(:Y, "Array of y-coordinates"),
    KwArgInfo((:ls, :linestyle), "Line style", :solid, values=_line_style_list),
    KwArgInfo(:dash, "Dash pattern", Float64[]),
    KwArgInfo((:linecolor, :lc, :color), "Line linecolor", :default),
    KwArgInfo((:lw, :lineweight), "Line weight", 0.5, cond=:(lw>0)),
    KwArgInfo(:marker, "Marker shape", :none,  values=_marker_list),
    KwArgInfo((:markersize, :ms), "Marker size", 2.5, cond=:(markersize>0)),
    KwArgInfo((:markercolor, :mc), "Marker color", :white),
    KwArgInfo((:mscolor, :markerstrokecolor, :msc), "Marker stroke color", :default),
    KwArgInfo(:label, "Data series label in legend", ""),
    KwArgInfo(:tag, "Data series tag over line", ""),
    KwArgInfo(:tagpos, "Tag position", 0.5),
    KwArgInfo(:tagloc, "Tag location", :top, values=[:bottom, :top, :left, :right]),
    KwArgInfo((:tagalong, :tagalign), "Sets that the tag will be aligned with the data series", false),
    KwArgInfo(:x, "x coordinate for a vertical line", nothing),
    KwArgInfo(:y, "y coordinate for a horizontal line", nothing),
    KwArgInfo(:order, "Order fo drawing", nothing),
    ArgCond(:(length(X)==length(Y)), "Length of X and Y arrays must be equal"),
]
@doc docstring(LineSeries_params) LineSeries

mutable struct LineSeries<:DataSeries
    X     ::AbstractArray
    Y     ::AbstractArray
    x     ::Union{Float64,Nothing}
    y     ::Union{Float64,Nothing}
    ls    ::Symbol
    lw    ::Float64
    linecolor::Union{Symbol,Tuple}
    marker::Symbol
    markersize::Float64
    markercolor::Union{Symbol,Tuple}
    mscolor::Union{Symbol,Tuple}
    label ::AbstractString
    tag::AbstractString
    tagpos::Float64
    tagloc::Symbol
    tagalong::Bool
    dash  ::Vector{Float64}
    order::Union{Int,Nothing}

    function LineSeries(X::AbstractArray, Y::AbstractArray; args...)

        args = checkargs([X, Y], args, LineSeries_params, aliens=false)

        if args.linecolor!==:default
            linecolor = get_color(args.linecolor)
            markercolor = get_color(args.markercolor, linecolor)
            mscolor = get_color(args.mscolor, linecolor)
        else
            linecolor = args.linecolor
            markercolor = args.markercolor
            mscolor = args.mscolor
        end

        n = min(length(X), length(Y))

        lw = args.lw
        ls = args.ls
        dash = args.dash

        if length(dash)==0
            if ls==:dash
                dash = [4.0, 2.4]*lw
            elseif ls==:dashdot
                dash = [2.0, 1.0, 2.0, 1.0]*lw
            elseif ls==:dot
                dash = [1.0, 1.0]*lw
            end
        else
            ls = :dash
        end

        this             = new(X[1:n], Y[1:n])
        this.ls          = ls
        this.lw          = lw
        this.linecolor   = linecolor
        this.marker      = args.marker
        this.markersize  = args.markersize
        this.markercolor = markercolor
        this.mscolor     = mscolor
        this.label       = args.label
        this.tag         = args.tag
        this.tagloc      = args.tagloc
        this.tagpos      = args.tagpos
        this.tagalong    = args.tagalong
        this.dash        = dash
        this.x           = args.x
        this.y           = args.y
        this.order       = args.order
        return this
    end
end


function LineSeries(; args...)
    return LineSeries(Float64[], Float64[]; args...)
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
    new_path(cc)

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