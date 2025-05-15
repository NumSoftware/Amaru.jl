# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct Colorbar<:FigureComponent
    location::Symbol
    colormap::Colormap
    axis::Axis
    thickness::Float64  # only bar thickness 
    scale::Float64
    inner_sep::Float64
    box::Vector{Float64}
    width::Float64
    height::Float64
    args::NamedTuple

    function Colorbar(; args...)
        args = checkargs(args, 
            [
                KwArgInfo( :location, "Colorbar location", :right, values=(:none, :right, :bottom) ),
                KwArgInfo( :colormap, "Colormap", :coolwarm),
                KwArgInfo( :limits, "Colorbar limit values", [0.0,0.0], length=2 ),
                KwArgInfo( :label, "Colorbar label", "", type=AbstractString ),
                KwArgInfo( :font_size, "Font size", 9.0, cond=:(font_size>0)),
                KwArgInfo( :font, "Font name", "NewComputerModern", type=AbstractString),
                KwArgInfo( :ticks, "Colorbar tick values", Float64[], type=AbstractArray ),
                KwArgInfo( :tick_labels, "Colorbar tick labels", String[], type=AbstractArray ),
                KwArgInfo( :tick_length, "Colorbar tick length", 3 ),
                KwArgInfo( :bins, "Number of bins", 6 ),
                KwArgInfo( :inner_sep, "Colorbar inner separation", 3 ),
                KwArgInfo( :scale, "Size scale", 1.0 ),
            ],
            aliens=false,
        )

        colormap = args.colormap isa Symbol ? Colormap(args.colormap) : args.colormap
        this = new(args.location, colormap)

        if args.location != :none
            direction = args.location == :right ? :vertical : :horizontal
            this.axis =  Axis(;
                direction   = direction,
                location    = args.location,
                limits      = args.limits,
                label       = args.label,
                font_size   = args.font_size,
                font        = args.font,
                bins        = args.bins,
                # ticks       = c.args.colorbarticks,
                # tick_labels = c.args.colorbarticklabels,
                # mult        = c.args.colorbarmult,
            )
        end

        this.scale = args.scale
        this.args = args
        return this
    end
end


function configure!(fig::Figure, cb::Colorbar)
    if cb.location!==:none
        configure!(cb.axis)
        # cb.thickness = 0.035*max(fig.width, fig.height)
        cb.thickness = 1.33*cb.args.font_size
        cb.inner_sep = cb.thickness
        if cb.location==:right
            cb.height = cb.scale*(fig.height - 2*fig.outerpad)
            cb.axis.height = cb.height
            cb.width = cb.inner_sep + cb.thickness + cb.axis.tick_length + cb.axis.width 
        elseif cb.location==:bottom
            cb.width = cb.scale*(fig.width - 2*fig.outerpad)
            cb.axis.width = cb.width
            cb.height = cb.inner_sep + cb.thickness + cb.axis.tick_length + cb.axis.height
        end
    end
end


function draw!(fig::Figure, cc::CairoContext, cb::Colorbar)
    set_matrix(cc, CairoMatrix([1, 0, 0, 1, 0, 0]...))

    if cb.location==:right
        # Axis
        fmin, fmax = cb.axis.limits
        x = fig.canvas.box[3]+cb.inner_sep+cb.thickness+cb.axis.tick_length
        h = cb.height

        y = fig.height/2 - h/2
        move_to(cc, x, y)
        draw!(cc, cb.axis)
        
        # Colorbar
        x = fig.canvas.box[3] + cb.inner_sep
        w = cb.thickness
        y = fig.height/2 + h/2

        pat = pattern_create_linear(0.0, y,  0.0, y-h)
        nstops = length(cb.colormap.stops)
        for i in 1:nstops
            stop = cb.colormap.stops[i]
            color = cb.colormap.colors[i]
            stop = round((stop-fmin)/(fmax-fmin), digits=8)
            pattern_add_color_stop_rgb(pat, stop, color...)
        end

        set_source(cc, pat)
        rectangle(cc, x, y, w, -h)
        fill(cc)
    elseif cb.location==:bottom
        # Axis
        fmin, fmax = cb.axis.limits
        w = cb.width
        x = fig.width/2 - w/2
        y = fig.canvas.box[4]+cb.inner_sep+cb.thickness+cb.axis.tick_length
        move_to(cc, x, y)
        draw!(cc, cb.axis)
        
        # Colorbar
        x = fig.width/2 - w/2
        y = fig.canvas.box[4] + cb.inner_sep
        h = cb.thickness

        pat = pattern_create_linear(x, 0.0,  x+w, 0.0)
        nstops = length(cb.colormap.stops)
        for i in 1:nstops
            stop = cb.colormap.stops[i]
            color = cb.colormap.colors[i]
            stop = round((stop-fmin)/(fmax-fmin), digits=8)
            pattern_add_color_stop_rgb(pat, stop, color...)
        end

        set_source(cc, pat)
        rectangle(cc, x, y, w, h)
        fill(cc)
    end
end