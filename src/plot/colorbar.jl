# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct Colorbar<:FigureComponent
    location::Symbol
    colormap::Colormap
    axis::Axis
    thickness::Float64  # only bar thickness 
    scale::Float64
    innersep::Float64
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
                KwArgInfo( :fontsize, "Font size", 9.0, cond=:(fontsize>0)),
                KwArgInfo( :font, "Font name", "NewComputerModern", type=AbstractString),
                KwArgInfo( :ticks, "Colorbar tick values", Float64[], type=AbstractArray ),
                KwArgInfo( :ticklabels, "Colorbar tick labels", String[], type=AbstractArray ),
                KwArgInfo( :ticklength, "Colorbar tick length", 3 ),
                KwArgInfo( :bins, "Number of bins", 6 ),
                KwArgInfo( :innersep, "Colorbar inner separation", 3 ),
                KwArgInfo( :scale, "Colorbar inner separation", 1.0 ),
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
                fontsize    = args.fontsize,
                font        = args.font,
                # ticks       = c.args.colorbarticks,
                # ticklabels  = c.args.colorbarticklabels,
                # mult        = c.args.colorbarmult,
            )
        end

        this.scale = args.scale
        this.args = args
        return this
    end
end


function configure!(c::Figure, cb::Colorbar)
    if cb.location!==:none
        configure!(cb.axis)
        cb.thickness = 0.035*max(c.width, c.height)
        cb.innersep = cb.thickness
        if cb.location==:right
            cb.height = cb.scale*(c.height - 2*c.outerpad)
            cb.axis.height = cb.height
            cb.width = cb.innersep + cb.thickness + cb.axis.ticklength + cb.axis.width 
        elseif cb.location==:bottom
            cb.width = cb.scale*(c.width - 2*c.outerpad)
            cb.axis.width = cb.width
            cb.height = cb.innersep + cb.thickness + cb.axis.ticklength + cb.axis.height
        end
    end
end


function draw!(fig::Figure, cc::CairoContext, cb::Colorbar)
    set_matrix(cc, CairoMatrix([1, 0, 0, 1, 0, 0]...))

    if cb.location==:right
        # Axis
        fmin, fmax = cb.axis.limits
        x = fig.canvas.box[3]+cb.innersep+cb.thickness+cb.axis.ticklength
        h = cb.height

        y = fig.height/2 - h/2
        move_to(cc, x, y)
        draw!(cc, cb.axis)
        
        # Colorbar
        x = fig.canvas.box[3] + cb.innersep
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
        y = fig.canvas.box[4]+cb.innersep+cb.thickness+cb.axis.ticklength
        move_to(cc, x, y)
        draw!(fig, cc, cb.axis)
        
        # Colorbar
        x = fig.width/2 - w/2
        y = fig.canvas.box[4] + cb.innersep
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