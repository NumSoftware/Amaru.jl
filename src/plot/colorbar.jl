# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct Colorbar<:ChartComponent
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
                ArgInfo( :location, "Colorbar location", :right, values=(:right, :bottom) ),
                ArgInfo( :colormap, "Colormap", :coolwarm),
                ArgInfo( :limits, "Colorbar limit values", [0.0,0.0], length=2 ),
                ArgInfo( :label, "Colorbar label", "", type=AbstractString ),
                ArgInfo( :fontsize, "Font size", 9.0, condition=:(fontsize>0)),
                ArgInfo( :font, "Font name", "NewComputerModern", type=AbstractString),
                ArgInfo( :ticks, "Colorbar tick values", Float64[], type=AbstractArray ),
                ArgInfo( :ticklabels, "Colorbar tick labels", String[], type=AbstractArray ),
                ArgInfo( :ticklength, "Colorbar tick length", 3 ),
                ArgInfo( :bins, "Number of bins", 6 ),
                ArgInfo( :innersep, "Colorbar inner separation", 3 ),
                ArgInfo( :scale, "Colorbar inner separation", 1.0 ),
            ],
            aliens=false,
        )

        colormap = args.colormap isa Symbol ? Colormap(args.colormap) : args.colormap
        this = new(args.location, colormap)

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

        this.scale = args.scale
        this.args = args
        return this
    end
end


function configure!(c::AbstractChart, cb::Colorbar)
    configure!(c, cb.axis)
    cb.thickness = 0.05*minimum(c.figsize)
    cb.innersep = cb.thickness
    if cb.location==:right
        cb.height = cb.scale*(c.figsize[2] - 2*c.outerpad)
        cb.axis.height = cb.height
        cb.width = cb.innersep + cb.thickness + cb.axis.ticklength + cb.axis.width 
    else
        cb.width = cb.scale*(c.figsize[1] - 2*c.outerpad)
        cb.axis.width = cb.width
        cb.height = cb.innersep + cb.thickness + cb.axis.ticklength + cb.axis.height
    end
end


function draw!(c::AbstractChart, cc::CairoContext, cb::Colorbar)
    fmin, fmax = cb.axis.limits
    set_matrix(cc, CairoMatrix([1, 0, 0, 1, 0, 0]...))

    if cb.location==:right
        # Axis
        x = c.canvas.box[3]+cb.innersep+cb.thickness+cb.axis.ticklength
        h = cb.height

        y = c.figsize[2]/2 - h/2
        move_to(cc, x, y)
        draw!(c, cc, cb.axis)
        
        # Colorbar
        x = c.canvas.box[3] + cb.innersep
        w = cb.thickness
        y = c.figsize[2]/2 + h/2

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
    else
        # Axis
        w = cb.width
        x = c.figsize[1]/2 - w/2
        y = c.canvas.box[4]+cb.innersep+cb.thickness+cb.axis.ticklength
        move_to(cc, x, y)
        draw!(c, cc, cb.axis)
        
        # Colorbar
        x = c.figsize[1]/2 - w/2
        y = c.canvas.box[4] + cb.innersep
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