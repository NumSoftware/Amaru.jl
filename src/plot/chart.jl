# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

abstract type AbstractChart end
abstract type ChartComponent end
abstract type DataSeriesPlot end

const _legend_positions=[
    :right,
    :left,
    :top,
    :bottom,
    :topright,
    :topleft,
    :bottomright,
    :bottomleft,
    :outerright,
    :outerleft,
    :outertopright,
    :outertopleft,
    :outerbottomright,
    :outerbottomleft
]


mutable struct Chart<:AbstractChart
    figsize::Union{Array, Tuple}
    xaxis::Union{ChartComponent, Nothing}
    yaxis::Union{ChartComponent, Nothing}
    canvas::Union{ChartComponent, Nothing}
    legend::Union{ChartComponent, Nothing}
    colorbar::Union{ChartComponent, Nothing}
    dataseries::Array

    outerpad::Float64
    icolor::Int
    args::NamedTuple

    function Chart(; args...)

        args = checkargs( args, 
            ArgInfo( :figsize, "Chart drawing size in dpi", default=(300,200), length=2),
            ArgInfo( :fontsize, "Font size", default=7.0, condition=:(fontsize>0)),
            ArgInfo( :xlimits, "x-axis limit values", default=[0.0,0.0], length=2 ),
            ArgInfo( :ylimits, "y-axis limit values", default=[0.0,0.0], length=2 ),
            ArgInfo( :xmult, "x-axis values multiplier", default=1.0 ),
            ArgInfo( :ymult, "y-axis values multiplier", default=1.0 ),
            ArgInfo( :xbins, "Number of bins in the x axis", default=7 ),
            ArgInfo( :ybins, "Number of bins in the y axis", default=6 ),
            ArgInfo( :xlabel, "Label for the x axis", default="", type=AbstractString ),
            ArgInfo( :ylabel, "Label for the y axis", default="", type=AbstractString ),
            ArgInfo( :xticks, "x-axis tick values", default=Float64[], type=AbstractArray ),
            ArgInfo( :yticks, "y-axis tick values", default=Float64[], type=AbstractArray ),
            ArgInfo( :xticklabels, "x-axis tick labels", default=String[], type=AbstractArray ),
            ArgInfo( :yticklabels, "y-axis tick labels", default=String[], type=AbstractArray ),
            ArgInfo( (:legendloc, :legend), "Legend location", default=:topright, values=_legend_positions ),
            ArgInfo( :legendfontsize, "Legend font size", default=7.0, condition=:(legendfontsize>0)),
            ArgInfo( (:colorbarloc, :colorbar), "Colorbar location", default=:right, values=(:right, :bottom) ),
            ArgInfo( (:colorbarscale, :cbscale), "Colorbar scale", default=0.9, condition=:(colorbarscale>0) ),
            ArgInfo( (:colorbarlabel, :cblabel, :colorbartitle), "Colorbar label", default="" ),
            ArgInfo( (:colorbarlimits, :cblimits), "Colorbar limits", default=Float64[0.0,0.0], length=2 ),
            ArgInfo( (:colorbarfontsize, :cbfontsize), "Colorbar font size", default=7.0, condition=:(colorbarfontsize>0)),
        )

        this = new()
        this.figsize = args.figsize
        this.xaxis = nothing
        this.yaxis = nothing
        this.legend = nothing
        this.colorbar = nothing
        this.dataseries = []
        this.icolor = 1
        this.args = args

        return this
    end
end


function configure!(c::Chart)

    width, height = c.figsize
    c.outerpad = 0.01*minimum(c.figsize)

    c.xaxis = Axis(; 
        direction   = :horizontal,
        limits      = c.args.xlimits,
        label       = c.args.xlabel,
        fontsize    = c.args.fontsize,
        ticks       = c.args.xticks,
        ticklabels  = c.args.xticklabels,
        mult        = c.args.xmult,
    )
    
    c.yaxis = Axis(; 
        direction   = :vertical,
        limits      = c.args.ylimits,
        label       = c.args.ylabel,
        fontsize    = c.args.fontsize,
        ticks       = c.args.yticks,
        ticklabels  = c.args.yticklabels,
        mult        = c.args.ymult,
    )

    configure!(c, c.xaxis)
    configure!(c, c.yaxis)

    # set width and height of axes
    c.xaxis.width  = width - c.yaxis.width - 2*c.outerpad
    c.yaxis.height = height - c.xaxis.height - 2*c.outerpad

    # set width and height of c
    c.canvas = Canvas()
    c.canvas.width = width - c.yaxis.width - 2*c.outerpad
    c.canvas.height = height - c.xaxis.height - 2*c.outerpad
    c.canvas.box = [ c.outerpad + c.yaxis.width, c.outerpad, width-c.outerpad, height - c.xaxis.height-c.outerpad ]


    c.canvas.limits = [ c.xaxis.limits[1], c.yaxis.limits[1], c.xaxis.limits[2], c.yaxis.limits[2] ]

    for p in c.dataseries
        configure!(c, p)
    end

    c.legend = Legend(; location=c.args.legendloc, fontsize=c.args.legendfontsize)
    configure!(c, c.legend)

end

function draw!(c::Chart, cc::CairoContext)
    set_source_rgb(cc, 1, 1, 1)  # White color
    paint(cc)

    x = c.outerpad+c.yaxis.width
    y = c.outerpad+c.yaxis.height
    move_to(cc, x, y)
    draw!(c, cc, c.xaxis)
    
    x = c.outerpad
    y = c.outerpad
    move_to(cc, x, y)
    draw!(c, cc, c.yaxis)

    # draw canvas grid
    draw!(c, cc, c.canvas)

    # draw plots
    x, y = c.canvas.box[1:2]
    w, h = c.canvas.box[3:4] - c.canvas.box[1:2]
    rectangle(cc, x, y, w, h)
    clip(cc)

    for p in c.dataseries
        draw!(c, cc, p)
    end
    reset_clip(cc)

    draw!(c, cc, c.legend)
end


function addplot!(c::Chart, P::DataSeriesPlot...)
    for p in P
        push!(c.dataseries, p)
    end
end


function save(chart::Chart, filename::String)
    width, height = chart.figsize
    surf = CairoPDFSurface(filename, width, height)
    cc = CairoContext(surf)
    
    configure!(chart)
    draw!(chart, cc)
    
    finish(surf)
end
