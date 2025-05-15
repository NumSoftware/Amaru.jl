# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct Axis<:FigureComponent
    direction::Symbol
    location::Symbol
    limits::Vector{Float64}
    label::AbstractString
    font::String
    font_size::Float64
    ticks::AbstractArray
    tick_labels::AbstractArray
    tick_length::Union{Float64,Nothing}
    n_bins::Int
    mult::Float64
    inner_sep::Float64
    width::Float64
    height::Float64

    function Axis(; args...)
        args = checkargs(args, func_params(Axis), aliens=false)

        if args.location==:none
            if args.direction==:horizontal
                location = :bottom
            else
                location = :left
            end
        else
            location = args.location
        end

        if length(args.tick_labels)>0
            length(args.ticks)!=length(args.tick_labels) && throw(AmaruException("Axis: the length of labels must match the number of ticks"))
        end

        limits = collect(args.limits)

        return new(args.direction, location, limits, args.label, args.font, args.font_size, args.ticks, args.tick_labels, args.tick_length, args.bins, args.mult, 3)
    end
end

func_params(::Type{Axis}) = [
    FunInfo( :Axis, "Creates an instance of an `Axis`."),
    KwArgInfo( :direction, "Axis direction", :horizontal, values=(:horizontal, :vertical) ),
    KwArgInfo( :location, "Axis location", :none, values=(:none, :left, :right, :top, :bottom) ),
    KwArgInfo( :limits, "Axis limit values", [0.0,0.0], length=2 ),
    KwArgInfo( :label, "Axis label", "", type=AbstractString ),
    KwArgInfo( :font, "Font name", "NewComputerModern", type=AbstractString),
    KwArgInfo( :font_size, "Font size", 7.0, cond=:(font_size>0)),
    KwArgInfo( :ticks, "Axis tick values", Float64[], type=AbstractArray ),
    KwArgInfo( :tick_labels, "Axis tick labels", String[], type=AbstractArray ),
    KwArgInfo( :tick_length, "Axis tick length", 3 ),
    KwArgInfo( :bins, "Hint for the number of bins", 6 ),
    KwArgInfo( :mult, "Axis values multiplier", 1.0 ),
    KwArgInfo( :inner_sep, "Axis inner pad", 3 ),
]
@doc docstring(func_params(Axis)) Axis()


function get_bin_length(vinf, vsup, n)
    dx = (vsup-vinf)/n*sign(vsup-vinf)
    ex = floor(log10(abs(dx))) # exponent 
    mant = dx/10^ex

    if mant<=1.4
        mant=1
    elseif mant<=1.8
        mant=1.5
    elseif mant<=2.25
        mant=2.0
    elseif mant<=3.0
        mant=2.5
    elseif mant<=4.5
        mant=4.0
    elseif mant<=7.5
        mant=5.0
    else
        mant=10.0
    end

    return round(mant*10^ex, sigdigits=2)*sign(vsup-vinf)
end


function make_ticklabels(ticks)
    # find mantissas and exponents
    M = Float64[]  # mantissas
    E = Int[]      # exponents

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
    max_digits = min(max_digits, 5)

    labels = []
    for (m, ex) in zip(M, E)
        if ex==0
            label = Printf.format(Printf.Format("%."*string(max_digits)*"f"), m)
        else
            # label = Printf.format(Printf.Format("%."*string(max_digits)*"f"), m)
            label = string(m)
            label *= "\\times 10^{$ex}"
        end
        push!(labels, latexstring(label))
    end

    return labels
end


function configure!(ax::Axis)

    # configure ticks
    if length(ax.ticks)==0
        len  = diff(ax.limits)[1]
        vinf, vsup = ax.limits

        if len==0
            vinf -= eps()*sign(vsup-vinf)
            vsup += eps()*sign(vsup-vinf)
        end
        dv = get_bin_length(vinf, vsup, ax.n_bins)

        # roundup first tick
        m = mod(vinf,dv)
        vinf = m==0 ? vinf : vinf - m + dv

        ax.ticks = round.(vinf:dv:vsup, digits=10)
    else # check ticks
        vmin, = minimum(ax.limits)
        vmax, = maximum(ax.limits)
        idxs = [ i for (i,tick) in enumerate(ax.ticks) if vmin<=tick<=vmax ]
        ax.ticks = ax.ticks[idxs]
        
        length(ax.tick_labels)!=0 && (ax.tick_labels = ax.tick_labels[idxs])
    end

    if length(ax.tick_labels)!=length(ax.ticks)
        ax.tick_labels = make_ticklabels(ax.ticks)
    end

    # configure latex strings
    if !isa(ax.label, LaTeXString) && contains(ax.label, "\$")
        ax.label = LaTeXString(ax.label)
    end

    ax.tick_labels = [ !isa(label, LaTeXString) && contains(label, "\$") ? LaTeXString(label) : label for label in ax.tick_labels  ]
    ax.n_bins = length(ax.ticks) - 1

    # configure size (axis dimensions do do not include tick lengths)
    ax.tick_length = 0.4*ax.font_size
    
    if ax.direction == :horizontal
        ax.inner_sep = 0.4*ax.font_size
        tk_lbs_height = maximum( getsize(lbl, ax.font_size)[2] for lbl in ax.tick_labels )
        label_height = getsize(ax.label, ax.font_size)[2]
        ax.height = label_height + ax.inner_sep + tk_lbs_height + ax.tick_length
    else
        tk_lbs_width = maximum( getsize(lbl, ax.font_size)[1] for lbl in ax.tick_labels )
        label_height = getsize(ax.label, ax.font_size)[2]
        ax.width = label_height + ax.inner_sep + tk_lbs_width + ax.tick_length
    end
    
    return ax
end



function draw!(cc::CairoContext, ax::Axis)

    Cairo.save(cc)

    x0, y0 = get_current_point(cc)
    
    font = get_font(ax.font)
    select_font_face(cc, font, Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL )
    set_font_size(cc, ax.font_size)
    set_matrix(cc, CairoMatrix([1, 0, 0, 1, 0, 0]...))
    
    set_source_rgb(cc, 0, 0, 0) # black
    set_line_width(cc, 0.05*ax.font_size)
    
    if ax.direction==:horizontal
        tk_lbs_height = maximum( getsize(lbl, ax.font_size)[2] for lbl in ax.tick_labels )
        label_height = getsize(ax.label, ax.font_size)[2]
        xmin, xmax = ax.limits

        # draw tick labels
        for (x,label) in zip(ax.ticks, ax.tick_labels)
            min(xmin, xmax) <= x <=max(xmin, xmax) || continue
            x1 = x0 + ax.width/(xmax-xmin)*(x-xmin)

            move_to(cc, x1, y0); rel_line_to(cc, 0, -ax.tick_length); stroke(cc)
            draw_text(cc, x1, y0+ax.tick_length+tk_lbs_height/2, label, halign="center", valign="center")
        end

        x = x0 + ax.width/2
        y = y0 + ax.height - label_height/2
        # y = y0 + ax.tick_length + ax.inner_sep + label_height
        draw_text(cc, x, y, ax.label, halign="center", valign="center", angle=0)

    else # :vertical ax
        tk_lbs_width = maximum( getsize(lbl, ax.font_size)[1] for lbl in ax.tick_labels )
        label_height = getsize(ax.label, ax.font_size)[2]
        ymin, ymax = ax.limits

        if ax.location==:left
            halign="right"
            tick_length = ax.tick_length
            x1 = x0 + ax.width
        else
            halign="left"
            tick_length = -ax.tick_length
            x1 = x0
        end

        # draw tick labels
        for (y,label) in zip(ax.ticks, ax.tick_labels)
            min(ymin, ymax) <= y <=max(ymin, ymax) || continue
            y1 = y0 + ax.height/(ymax-ymin)*(ymax-y)
            
            move_to(cc, x1, y1); rel_line_to(cc, tick_length, 0); stroke(cc)

            draw_text(cc, x1-tick_length, y1+ax.font_size/2, label, halign=halign, valign="bottom")
        end
        
        # draw label
        if ax.location==:left
            # x = x0 + label_height/2
            x = x0 + ax.width - ax.tick_length - tk_lbs_width - ax.inner_sep - label_height/2
            # x = x0 + ax.width - tick_length - tk_lbs_width - tick_length 
        else
            x = x0 + ax.tick_length + tk_lbs_width + ax.inner_sep + label_height/2
        end
        y = y0 + ax.height/2

        draw_text(cc, x, y, ax.label, halign="center", valign="center", angle=90)
    end

    Cairo.restore(cc)
end