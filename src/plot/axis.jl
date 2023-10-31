# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct Axis<:ChartComponent
    direction::Symbol
    location::Symbol
    limits::Vector{Float64}
    label::AbstractString
    fontsize::Float64
    ticks::Array
    ticklabels::Array
    ticklength::Union{Float64,Nothing}
    nbins::Int
    mult::Float64
    innersep::Float64
    width::Float64
    height::Float64

    function Axis(; args...)
        args = checkargs(args,
            ArgInfo( :direction, "Axis direction", default=:horizontal, values=(:horizontal, :vertical) ),
            ArgInfo( :location, "Axis location", default=:none, values=(:none, :left, :right, :top, :bottom) ),
            ArgInfo( :limits, "Axis limit values", default=[0.0,0.0], length=2 ),
            ArgInfo( :label, "Axis label", default="", type=AbstractString ),
            ArgInfo( :fontsize, "Font size", default=7.0, condition=:(fontsize>0)),
            ArgInfo( :ticks, "Axis tick values", default=Float64[], type=AbstractArray ),
            ArgInfo( :ticklabels, "Axis tick labels", default=String[], type=AbstractArray ),
            ArgInfo( :ticklength, "Axis tick length", default=3 ),
            ArgInfo( :bins, "Number of bins", default=6 ),
            ArgInfo( :mult, "Axis values multiplier", default=1.0 ),
            ArgInfo( :innersep, "Axis inner pad", default=3 ),
        )

        if args.location==:none
            if args.direction==:horizontal
                location = :bottom
            else
                location = :left
            end
        else
            location = args.location
        end

        return new(args.direction, location, args.limits, args.label, args.fontsize, args.ticks, args.ticklabels, args.ticklength, args.bins, args.mult, 3)
    end
end


function roundbin(dx)
    ex = floor(log10(dx)) # exponent 
    mant = dx/10^ex

    if mant<=1.4
        mant=1
    elseif mant<=1.8
        mant=1.5
    elseif mant<=2.25
        mant=2.0
    elseif mant<=3.5
        mant=2.5
    elseif mant<=7.5
        mant=5.0
    else
        mant=10.0
    end

    return round(mant*10^ex, sigdigits=2)
end


function configure!(chart::AbstractChart, ax::Axis)
    ax.ticklength = 0.015*minimum(chart.figsize)
    ax.innersep   = 0.02*minimum(chart.figsize)

    # configure limits
    if ax.limits==[0.0,0.0]
        if length(ax.ticks)==0
            limits = [Inf, -Inf]
            for p in chart.dataseries
                p isa DataSeriesPlot || continue
                if ax.direction==:horizontal
                    limits[1] = min(limits[1], minimum(p.X))
                    limits[2] = max(limits[2], maximum(p.X))
                else
                    limits[1] = min(limits[1], minimum(p.Y))
                    limits[2] = max(limits[2], maximum(p.Y))
                end
            end
            ax.limits = limits
        else
            limits = extrema(chart.args.xticks)
            ax.limits = collect(limits)
        end
    end

    # configure ticks
    if length(ax.ticks)==0
        len  = diff(ax.limits)[1]
        vmin, vmax = ax.limits

        if len==0
            vmin -= eps()
            vmax += eps()
        end

        dv = roundbin((vmax-vmin)/ax.nbins)

        # roundup first tick
        m = mod(vmin,dv) 
        vmin = m==0 ? vmin : round(vmin - m + dv, sigdigits=2)

        ax.ticks = round.(vmin:dv:vmax, digits=10)
        # ax.ticks = round.(vmin:dv:vmax, sigdigits=3)
    end

    if length(ax.ticklabels)!=length(ax.ticks)
        ax.ticklabels = make_ticklabels(ax.ticks)
    end

    # configure latex strings
    if !isa(ax.label, LaTeXString) && contains(ax.label, "\$")
        ax.label = LaTeXString(ax.label)
    end

    ax.ticklabels = [ !isa(label, LaTeXString) && contains(label, "\$") ? LaTeXString(label) : label for label in ax.ticklabels  ]

    ax.nbins = length(ax.ticks) - 1

    # configure size
    if ax.direction == :horizontal
        tk_lbs_height = maximum( getsize(lbl, ax.fontsize)[2] for lbl in ax.ticklabels )
        label_height = getsize(ax.label, ax.fontsize)[2]
        ax.height = label_height + ax.innersep + tk_lbs_height + ax.innersep
    else
        tk_lbs_width = maximum( getsize(lbl, ax.fontsize)[1] for lbl in ax.ticklabels )
        label_height = getsize(ax.label, ax.fontsize)[2]
        ax.width = label_height + ax.innersep + tk_lbs_width + ax.innersep
    end
    
    return ax
end

function make_ticklabels(ticks)
    # find mantissas and exponents
    M = Float64[]
    E = Int[]
    for x in ticks
        if abs(x)<1e-10
            x = 0.0
            expo = 1
        else
            expo = floor(Int, log10(abs(x)))
        end

        if -5<expo<5
            mantissa = x
            expo = 0
        else
            mantissa = round(x/10.0^expo, sigdigits=3)
        end
        push!(M, mantissa)
        push!(E, expo)
    end

    # find max decimal digits in mantissas
    max_digits = 0
    for m in M
        digits = 0
        fraction = round(abs(m-trunc(m)), sigdigits=10)
        while fraction > 0
            digits += 1
            fraction *= 10
            fraction = round(abs(fraction-trunc(fraction)), sigdigits=10)
        end
        max_digits = max(max_digits, digits)
    end

    labels = []
    for (m, ex) in zip(M, E)
        label = Printf.format(Printf.Format("%."*string(max_digits)*"f"), m)
        ex!=0 && (label *= "\\times 10^{$ex}")
        push!(labels, latexstring(label))
    end

    return labels
end


function draw!(c::AbstractChart, cc::CairoContext, ax::Axis)
    x0, y0 = get_current_point(cc)

    set_font_size(cc, ax.fontsize)
    set_font_face(cc, "NewComputerModern $(ax.fontsize)") # for pango text
    set_matrix(cc, CairoMatrix([1, 0, 0, 1, 0, 0]...))
    
    set_source_rgb(cc, 0, 0, 0) # black
    set_line_width(cc, 0.5)
    
    if ax.direction==:horizontal
        tk_lbs_height = maximum( getsize(lbl, ax.fontsize)[2] for lbl in ax.ticklabels )
        label_height = getsize(ax.label, ax.fontsize)[2]
        xmin, xmax = ax.limits

        for (x,label) in zip(ax.ticks, ax.ticklabels)
            x1 = x0 + ax.width/(xmax-xmin)*(x-xmin)

            move_to(cc, x1, y0); line_to(cc, x1, y0-ax.ticklength); stroke(cc)
            if label isa LaTeXString
                textext(cc, x1, y0+ax.ticklength+tk_lbs_height/2, label, halign="center", valign="center")
            else
                text(cc, x1, y0+ax.ticklength+tk_lbs_height/2, label, halign="center", valign="center")
            end
        end

        x = x0 + ax.width/2
        y = y0 + ax.height - label_height/2

        if ax.label isa LaTeXString
            textext(cc, x, y, ax.label, halign="center", valign="center", angle=0)
        else
            text(cc, x, y, ax.label, halign="center", valign="center", angle=0)
        end

    else # vertical ax
        tk_lbs_width = maximum( getsize(lbl, ax.fontsize)[1] for lbl in ax.ticklabels )
        label_height = getsize(ax.label, ax.fontsize)[2]


        if ax.location==:left
            halign="right"
            ticklength = ax.ticklength
            x1 = x0 + ax.width
        else
            halign="left"
            ticklength = -ax.ticklength
            x1 = x0
        end

        ymin, ymax = ax.limits
        for (y,label) in zip(ax.ticks, ax.ticklabels)
            y1 = y0 + ax.height/(ymax-ymin)*(ymax-y)
            
            move_to(cc, x1, y1); line_to(cc, x1+ticklength, y1); stroke(cc)

            if label isa LaTeXString
                textext(cc, x1-ticklength, y1, label, halign=halign, valign="center")
            else
                text(cc, x1-ticklength, y1, label, halign=halign, valign="center")
            end
        end
        
        if ax.location==:left
            x = x0 + label_height/2
        else
            x = x0 + tk_lbs_width + ax.ticklength + label_height/2
        end
        y = y0 + ax.height/2

        if ax.label isa LaTeXString
            textext(cc, x, y, ax.label, halign="center", valign="center", angle=90)
        else
            text(cc, x, y, ax.label, halign="center", valign="center", angle=90)
        end
    end

end