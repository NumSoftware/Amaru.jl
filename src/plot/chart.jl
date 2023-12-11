# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

abstract type AbstractChart end
abstract type ChartComponent end
abstract type DataSeriesPlot end

_available_formats=[
    ".pdf",
    ".png",
    ".svg",
]

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
    width::Float64
    height::Float64
    xaxis::Union{ChartComponent, Nothing}
    yaxis::Union{ChartComponent, Nothing}
    canvas::Union{ChartComponent, Nothing}
    legend::Union{ChartComponent, Nothing}
    colorbar::Union{ChartComponent, Nothing}
    dataseries::Array

    outerpad::Float64
    toppad::Float64
    rightpad::Float64
    icolor::Int
    iorder::Int
    args::NamedTuple

    function Chart(; args...)
        args = checkargs(args, func_params(Chart), aliens=false)
            
        this = new()
        this.width, this.height = args.size
        this.xaxis = nothing
        this.yaxis = nothing
        this.legend = nothing
        this.colorbar = nothing
        this.dataseries = []
        this.icolor = 1
        this.iorder = 1
        this.args = args

        return this
    end
end


func_params(::Type{Chart}) = [
    FunInfo( :Chart, "Creates a customizable `Chart` instance.", ()),
    ArgInfo( (:size, :figsize), "Chart drawing size in dpi", (220,150), length=2),
    ArgInfo( :font, "Font name", "NewComputerModern", type=AbstractString),
    ArgInfo( :fontsize, "Font size", 7.0, condition=:(fontsize>0)),
    ArgInfo( :xlimits, "x-axis limit values", [0.0,0.0], length=2 ),
    ArgInfo( :ylimits, "y-axis limit values", [0.0,0.0], length=2 ),
    ArgInfo( :aspectratio, "aspect ratio", :auto, values=(:auto, :equal) ),
    ArgInfo( :xmult, "x-axis values multiplier", 1.0 ),
    ArgInfo( :ymult, "y-axis values multiplier", 1.0 ),
    ArgInfo( :xbins, "Number of bins in the x axis", 7 ),
    ArgInfo( :ybins, "Number of bins in the y axis", 6 ),
    ArgInfo( :xlabel, "Label for the x axis", L"$x$", type=AbstractString ),
    ArgInfo( :ylabel, "Label for the y axis", L"$y$", type=AbstractString ),
    ArgInfo( :xticks, "x-axis tick values", Float64[], type=AbstractArray ),
    ArgInfo( :yticks, "y-axis tick values", Float64[], type=AbstractArray ),
    ArgInfo( :xticklabels, "x-axis tick labels", String[], type=AbstractArray ),
    ArgInfo( :yticklabels, "y-axis tick labels", String[], type=AbstractArray ),
    ArgInfo( (:legendloc, :legend), "Legend location", :topright, values=_legend_positions ),
    ArgInfo( :legendfontsize, "Legend font size", :fontsize, condition=:(legendfontsize>0)),
    ArgInfo( (:colorbarloc, :colorbar), "Colorbar location", :right, values=(:right, :bottom) ),
    ArgInfo( (:colorbarscale, :cbscale), "Colorbar scale", 0.9, condition=:(colorbarscale>0) ),
    ArgInfo( (:colorbarlabel, :cblabel, :colorbartitle), "Colorbar label", "" ),
    ArgInfo( (:colorbarlimits, :cblimits), "Colorbar limits", Float64[0.0,0.0], length=2 ),
    ArgInfo( (:colorbarfontsize, :cbfontsize), "Colorbar font size", 7.0, condition=:(colorbarfontsize>0)),
]
@doc make_doc(Chart) Chart()


function configure!(c::Chart)

    # width, height = c.figsize
    c.outerpad = 0.01*min(c.width, c.height)
    c.toppad = c.outerpad
    c.rightpad = c.outerpad

    c.xaxis = Axis(; 
        direction  = :horizontal,
        limits     = c.args.xlimits,
        label      = c.args.xlabel,
        font       = c.args.font,
        fontsize   = c.args.fontsize,
        ticks      = c.args.xticks,
        ticklabels = c.args.xticklabels,
        mult       = c.args.xmult,
    )
    
    c.yaxis = Axis(; 
        direction  = :vertical,
        limits     = c.args.ylimits,
        label      = c.args.ylabel,
        font       = c.args.font,
        fontsize   = c.args.fontsize,
        ticks      = c.args.yticks,
        ticklabels = c.args.yticklabels,
        mult       = c.args.ymult,
    )

    # configure axes, may change chart pads
    configure!(c, c.xaxis, c.yaxis)
    # configure!(c, c.xaxis)
    # configure!(c, c.yaxis)

    # set width and height of canvas
    c.canvas = Canvas()

    configure!(c, c.canvas)


    for p in c.dataseries
        configure!(c, p)
    end

    has_legend = any( ds.label!="" for ds in c.dataseries )
    if has_legend
        # fontsize = c.args.legendfontsize > c.args.fontsize ? c.args.fontsize || c.args.legendfontsize
        c.legend = Legend(; 
            location = c.args.legendloc,
            font     = c.args.font,
            fontsize = c.args.legendfontsize
        )
        configure!(c, c.legend)
    end
end

function draw!(c::Chart, cc::CairoContext)
    # draw canvas grid
    draw!(c, cc, c.canvas)

    # draw axes
    x = c.outerpad+c.yaxis.width
    y = c.toppad+c.yaxis.height
    move_to(cc, x, y)
    draw!(c, cc, c.xaxis)
    
    x = c.outerpad
    y = c.toppad
    move_to(cc, x, y)
    draw!(c, cc, c.yaxis)

    # draw plots
    x, y = c.canvas.box[1:2]
    w, h = c.canvas.box[3:4] - c.canvas.box[1:2]
    rectangle(cc, x, y, w, h)
    Cairo.clip(cc)

    sorted = sort(c.dataseries, by=x->x.order)
    for p in sorted
        draw!(c, cc, p)
    end
    reset_clip(cc)

    has_legend = any( ds.label!="" for ds in c.dataseries )
    if has_legend
        draw!(c, cc, c.legend)
    end
end


function addplot!(chart::Chart, P::DataSeriesPlot...)
    length(P)>0 || throw(AmaruException("No dataseries added"))
    for p in P
        if p.linecolor===:default # update colors
            p.linecolor = _colors_dict[_default_colors[chart.icolor]]
            chart.icolor = mod(chart.icolor, length(_default_colors)) + 1
        end

        if p.order===nothing
            p.order = chart.iorder
            chart.iorder += 1
        end

        push!(chart.dataseries, p)
    end

end

function addplot!(c::Chart, plots::Array{<:DataSeriesPlot,1})
    addplot!(c, plots...)
end


function save(chart::Chart, filename::String, copypath::String="")
    width, height = chart.width, chart.height

    fmt = splitext(filename)[end]
    if fmt==".pdf"
        surf = CairoPDFSurface(filename, width, height)
    elseif fmt==".svg"
        surf = CairoSVGSurface(filename, width, height)
    elseif fmt==".ps"
        surf = CairoPSSurface(filename, width, height)
    elseif fmt==".png"
        surf = CairoImageSurface(width, height, Cairo.FORMAT_ARGB32)
    else
        formats = join(_available_formats, ", ", " and ")
        throw(AmaruException("Cannot save image to format $fmt. Available formats are: $formats"))
    end

    cc = CairoContext(surf)    
    configure!(chart)

    if fmt==".png"
        set_source_rgb(cc, 1.0, 1.0, 1.0) # RGB values for white
        paint(cc)
    end

    draw!(chart, cc)
    
    if fmt==".png"
        write_to_png(surf, filename)
    else
        finish(surf)
    end

    copypath!="" && cp(filename, joinpath(dirname(copypath), basename(filename)), force=true)

    return nothing
end
