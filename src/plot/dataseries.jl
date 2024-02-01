# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

const _line_style_list = [:none, :solid, :dot, :dash, :dashdot]
const _marker_list = [:none, :circle, :square, :triangle, :utriangle, :cross, :xcross, :diamond, :pentagon, :hexagon, :star]


mutable struct LinePlot<:DataSeriesPlot
    X     ::Array
    Y     ::Array
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

    function LinePlot(X::AbstractArray, Y::AbstractArray; args...)

        args = checkargs(args, func_params(LinePlot), aliens=false)

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


function LinePlot(; args...)
    return LinePlot(Float64[], Float64[]; args...)
end

func_params(::Type{LinePlot}) = [
    FunInfo( :LinePlot, "Creates a customizable `LinePlot` instance.", "X, Y"),
    ArgInfo( (:ls, :linestyle), "Line style", :solid, values=_line_style_list ),
    ArgInfo( :dash, "Dash pattern", Float64[] ),
    ArgInfo( (:linecolor, :lc, :color), "Line linecolor", :default),
    ArgInfo( (:lw, :lineweight), "Line weight", 0.5, condition=:(lw>0) ),
    ArgInfo( :marker, "Marker shape", :none,  values=_marker_list ),
    ArgInfo( (:markersize, :ms), "Marker size", 2.5, condition=:(markersize>0) ),
    ArgInfo( (:markercolor, :mc), "Marker color", :white ),
    ArgInfo( (:mscolor, :markerstrokecolor, :msc), "Marker stroke color", :default ),
    ArgInfo( :label, "Data series label in legend", ""),
    ArgInfo( :tag, "Data series tag over line", ""),
    ArgInfo( :tagpos, "Tag position", 0.5),
    ArgInfo( :tagloc, "Tag location", :top, values=[:bottom, :top, :left, :right]),
    ArgInfo( (:tagalong, :tagalign), "Sets that the tag will be aligned with the data series", false),
    ArgInfo( :x, "x coordinate for a vertical line", nothing),
    ArgInfo( :y, "y coordinate for a horizontal line", nothing),
    ArgInfo( :order, "Order fo drawing", nothing),
]
@doc make_doc(func_params(LinePlot)) LinePlot()

function data2user(c::Chart, x, y)
    Xmin, Ymin, Xmax, Ymax = c.canvas.box
    xmin, ymin, xmax, ymax = c.canvas.limits
    xD = Xmin + (Xmax-Xmin)/(xmax-xmin)*(x-xmin)
    yD = Ymin + (Ymax-Ymin)/(ymax-ymin)*(ymax-y)
    return xD, yD
end

function configure!(chart::Chart, p::LinePlot)
    xmin, ymin, xmax, ymax = chart.canvas.limits
    if p.x !== nothing
        p.X = [ p.x, p.x ]
        p.Y = [ ymin, ymax ]
    elseif p.y !== nothing
        p.X = [ xmin, xmax ]
        p.Y = [ p.y, p.y ]
    end
    # if length(p.dash)==0
    #     if p.ls==:dash
    #         p.dash = [4.0, 2.4]*p.lw
    #     elseif p.ls==:dashdot
    #         p.dash = [2.0, 1.0, 2.0, 1.0]*p.lw
    #     elseif p.ls==:dot
    #         p.dash = [1.0, 1.0]*p.lw
    #     end
    # end
end

function draw!(chart::Chart, cc::CairoContext, p::LinePlot)

    p.markercolor = get_color(p.markercolor, p.linecolor)
    p.mscolor = get_color(p.mscolor, p.linecolor)

    set_matrix(cc, CairoMatrix([1, 0, 0, 1, 0, 0]...)) 
    set_source_rgb(cc, p.linecolor...)
    set_line_width(cc, p.lw)
    set_line_join(cc, Cairo.CAIRO_LINE_JOIN_ROUND)

    # Draw lines
    new_path(cc)
    n = length(p.X)
    X = p.X*chart.xaxis.mult
    Y = p.Y*chart.yaxis.mult
    
    if p.ls!==:none
        x1, y1 = data2user(chart, X[1], Y[1])

        if p.ls==:solid
            move_to(cc, x1, y1)
            for i in 2:n
                x, y = data2user(chart, X[i], Y[i])
                line_to(cc, x, y); 
            end
            stroke(cc)
        else # dashed
            len = sum(p.dash)
            offset = 0.0
            set_dash(cc, p.dash, offset)
            move_to(cc, x1, y1)
            for i in 2:n
                x, y = data2user(chart, X[i], Y[i])
                line_to(cc, x, y); 
                offset = mod(offset + norm((x1-x,y1-y)), len)
                set_dash(cc, p.dash, offset)
                x1, y1 = x, y
            end
            stroke(cc)
            set_dash(cc, Float64[])
        end
    end

    # Draw markers
    for (x,y) in zip(X, Y)
        x, y = data2user(chart, x, y)
        draw_marker(cc, x, y, p.marker, p.markersize, p.markercolor, p.mscolor)
    end

    # Draw tag
    if p.tag!=""
        len = 0.0
        L = [ len ] # lengths
        for i in 2:length(X)
            len += norm((X[i]-X[i-1], Y[i]-Y[i-1]))
            push!(L, len)
        end
        lpos = p.tagpos*len # length to position

        i = findfirst(z->z>lpos, L)
        i = min(i, length(L)-1)

        # location coordinates
        x  = X[i] + (lpos-L[i])/(L[i+1]-L[i])*(X[i+1]-X[i])
        y  = Y[i] + (lpos-L[i])/(L[i+1]-L[i])*(Y[i+1]-Y[i])
        
        # location coordinates in user units
        x, y = data2user(chart, x, y)
        x1, y1 = data2user(chart, X[i], Y[i])
        x2, y2 = data2user(chart, X[i+1], Y[i+1])
        α = -atand(y2-y1, x2-x1) # tilt

        # pads
        pad = chart.args.fontsize*0.3
        dx = pad*cosd(α)
        dy = pad*sind(α)

        # Default location "top"
        if p.tagloc==:top
            va = "bottom"
            if 0<α<= 90 || -180 <α<= -90 
                ha = "right" 
                dx, dy = -dy, -dx
            else
                ha = "left"
                dx, dy = dy, dx
            end
        else
            va = "top"
            if 0<α<=90 || -180<α<=-90
                ha = "left"
                dx, dy = dy, dx
            else
                ha = "right"
                dx, dy = -dy, dx
            end
        end

        if p.tagalong
            ha = "center"
            dx = 0.0
            dy = p.tagloc==:top ? -pad : 0.0
        else
            α = 0.0
        end

        set_font_size(cc, chart.args.fontsize*0.9)
        font = get_font(chart.args.font)
        select_font_face(cc, font, Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL )
        set_source_rgb(cc, 0, 0, 0)
        draw_text(cc, x+dx, y+dy, p.tag, halign=ha, valign=va, angle=α)
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