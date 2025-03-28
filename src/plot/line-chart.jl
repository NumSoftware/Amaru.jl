# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# abstract type Chart end
# abstract type FigureComponent end
# abstract type DataSeries end


Chart_params = [
    FunInfo( :Chart, "Creates a customizable `Chart` instance."),
    KwArgInfo( (:size, :figsize), "Chart drawing size in dpi", (220,150), length=2),
    KwArgInfo( :font, "Name of the font", "NewComputerModern", type=AbstractString),
    KwArgInfo( :fontsize, "Size of the font in dpi", 7.0, cond=:(fontsize>0)),
    KwArgInfo( (:xlimits, :xlims), "Limits of the x-axis", [0.0,0.0], length=2 ),
    KwArgInfo( (:ylimits, :ylims), "Limits of the y-axis", [0.0,0.0], length=2 ),
    KwArgInfo( :aspectratio, "Ratio of y-unit to x-unit", :auto, values=(:auto, :equal) ),
    KwArgInfo( :xmult, "Multiplier for x-axis values", 1.0 ),
    KwArgInfo( :ymult, "Multiplier for y-axis values", 1.0 ),
    KwArgInfo( :xbins, "Quantity of bins along the x-axis", 7 ),
    KwArgInfo( :ybins, "Quantity of bins along the y-axis", 6 ),
    KwArgInfo( :xlabel, "Text label for the x-axis", L"$x$", type=AbstractString ),
    KwArgInfo( :ylabel, "Text label for the y-axis", L"$y$", type=AbstractString ),
    KwArgInfo( :xticks, "Tick positions along the x-axis", Float64[], type=AbstractArray ),
    KwArgInfo( :yticks, "Tick positions along the y-axis", Float64[], type=AbstractArray ),
    KwArgInfo( :xticklabels, "Tick labels for each tick mark on the x-axis", String[], type=AbstractArray ),
    KwArgInfo( :yticklabels, "Tick labels for each tick mark on the y-axis", String[], type=AbstractArray ),
    KwArgInfo( (:legendloc, :legend), "Position of the legend on the chart", :topright, values=_legend_positions ),
    KwArgInfo( :legendfontsize, "Font size for the legend", :fontsize, cond=:(legendfontsize>0)),
    KwArgInfo( (:colorbarloc, :colorbar), "Placement of the colorbar", :right, values=(:right, :bottom) ),
    KwArgInfo( (:colorbarscale, :cbscale), "Scaling factor for the colorbar", 0.9, cond=:(colorbarscale>0) ),
    KwArgInfo( (:colorbarlabel, :cblabel), "Label for the colorbar", "" ),
    KwArgInfo( (:colorbarlimits, :cblimits), "Range of values in the colorbar", Float64[0.0,0.0], length=2 ),
    KwArgInfo( (:colorbarfontsize, :cbfontsize), "Font size for the colorbar labels", 7.0, cond=:(colorbarfontsize>0)),
]
@doc docstring(Chart_params) Chart


mutable struct Chart<:AbstractChart
    width::Float64
    height::Float64
    canvas::Union{Canvas, Nothing}
    xaxis::Union{Axis, Nothing}
    yaxis::Union{Axis, Nothing}
    legend::Union{Legend, Nothing}
    colorbar::Union{FigureComponent, Nothing}
    dataseries::Array{DataSeries,1}
    annotations::AbstractArray

    outerpad::Float64
    leftpad::Float64
    rightpad::Float64
    toppad::Float64
    bottompad::Float64
    icolor::Int
    iorder::Int
    args::NamedTuple

    function Chart(; args...)
        args = checkargs(args, Chart_params, aliens=false)
            
        this = new()
        this.width, this.height = args.size
        this.xaxis = nothing
        this.yaxis = nothing
        this.legend = nothing
        this.colorbar = nothing
        this.dataseries = []
        this.annotations = []
        this.icolor = 1
        this.iorder = 1
        this.args = args

        return this
    end
end


function addseries!(chart::Chart, series::DataSeries...)
    length(series)>0 || throw(AmaruException("No dataseries added"))
    for p in series
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


function addseries!(c::Chart, plots::Array{<:DataSeries,1})
    addseries!(c, plots...)
end


function addlegend!(c::Chart, legend::Legend) 
    c.legend = legend
end


function configure!(c::Chart)

    c.outerpad  = 0.01*min(c.width, c.height)
    c.leftpad   = c.outerpad
    c.rightpad  = c.outerpad
    c.toppad    = c.outerpad
    c.bottompad = c.outerpad

    # configure legend
    if c.legend===nothing
        # make legend if any dataseries has a label
        if any( ds.label!="" for ds in c.dataseries )
            # fontsize = c.args.legendfontsize > c.args.fontsize ? c.args.fontsize || c.args.legendfontsize
            c.legend = Legend(; 
                location = c.args.legendloc,
                font     = c.args.font,
                fontsize = c.args.legendfontsize,
                ncols    = 1
            )
            configure!(c, c.legend)
        end
    else
        configure!(c, c.legend)
    end

    # configure axes
    
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

    # configure the canvas
    configure!(c, c.canvas)

    for p in c.dataseries
        configure!(c, p)
    end

end


function configure!(c::Chart, canvas::Canvas)
    xmin, xmax = c.xaxis.limits
    ymin, ymax = c.yaxis.limits
    # canvas.limits = [ xmin, ymin, xmax, ymax ]

    if c.args.aspectratio==:equal
        # compute extra limits
        width = c.width - c.yaxis.width - c.leftpad - c.rightpad
        height = c.height - c.xaxis.height - c.toppad - c.rightpad
        r = min(width/(xmax-xmin), height/(ymax-ymin))
        dx = 0.5*(width/r - (xmax-xmin))
        dy = 0.5*(height/r - (ymax-ymin))
        
        # update limits
        c.xaxis.limits = [ xmin-dx, xmax+dx ]
        c.yaxis.limits = [ ymin-dy, ymax+dy ]

        # force recompute ticks
        c.xaxis.ticks = []
        c.yaxis.ticks = []

        # reconfigure axes
        configure!(c, c.xaxis, c.yaxis)
        xmin, xmax = c.xaxis.limits
        ymin, ymax = c.yaxis.limits

        # udpa
    end
    
    canvas.width = c.width - c.yaxis.width - c.leftpad - c.rightpad
    canvas.height = c.height - c.xaxis.height - c.toppad - c.bottompad
    canvas.box = [ c.leftpad + c.yaxis.width, c.toppad, c.width-c.rightpad, c.height - c.xaxis.height-c.bottompad ]
    canvas.limits = [ xmin, ymin, xmax, ymax ]
end




function configure!(chart::Chart, xax::Axis, yax::Axis)

    # check limits
    for ax in (xax, yax)
        if ax.limits==[0.0,0.0]
            if length(ax.ticks)==0
                limits = [Inf, -Inf]
                for p in chart.dataseries
                    p isa DataSeries || continue
                    if ax.direction==:horizontal
                        if p.x!==nothing
                            limits[1] = min(limits[1], p.x)
                            limits[2] = max(limits[2], p.x)
                        end
                        if length(p.X)>0
                            limits[1] = min(limits[1], minimum(p.X))
                            limits[2] = max(limits[2], maximum(p.X))
                        end
                    else
                        if p.y!==nothing
                            limits[1] = min(limits[1], p.y)
                            limits[2] = max(limits[2], p.y)
                        end
                        if length(p.Y)>0
                            limits[1] = min(limits[1], minimum(p.Y))
                            limits[2] = max(limits[2], maximum(p.Y))
                        end
                        
                    end
                end
                if limits[1]==limits[2]
                    limits = [-0.1, 0.1]
                end
            else
                limits = collect(extrema(chart.args.xticks))
            end

            # extend limits
            f = ax.direction == :horizontal ? 0.03 : 0.03*chart.width/chart.height
            dx = f*(limits[2]-limits[1])
            limits = [ limits[1]-dx, limits[2]+dx ]
            limits = [ limits[1]*ax.mult, limits[2]*ax.mult ]
            ax.limits = limits
        end
    end


    configure!(xax)
    configure!(yax)

    # set width and height of axes
    width, height = chart.width, chart.height
    xax.width  = width - yax.width - chart.leftpad - chart.rightpad
    yax.height = height - xax.height - chart.toppad - chart.bottompad

    # update chart.rightpad if required
    label_width = getsize(xax.ticklabels[end], xax.fontsize)[1] 
    xdist = xax.width*(xax.limits[2]-xax.ticks[end])/(xax.limits[2]-xax.limits[1]) # distance of the right most tick to the right side of axis
    if xdist-label_width/2 < 0
        chart.rightpad += label_width/2 - xdist
        xax.width = width - yax.width - chart.leftpad - chart.rightpad
    end

    # update chart.toppad if required
    label_height = getsize(yax.ticklabels[end], yax.fontsize)[2] 
    ydist = yax.height * (yax.limits[2] - yax.ticks[end])/(yax.limits[2]-yax.limits[1]) # distance of the upper most tick to the top side of axis
    if ydist-label_height/2 < 0
        chart.toppad += label_height/2 - ydist
        yax.height = height - xax.height - chart.toppad - chart.bottompad
    end

end


function configure!(chart::Chart, p::LineSeries)
    xmin, ymin, xmax, ymax = chart.canvas.limits
    if p.x !== nothing
        p.X = [ p.x, p.x ]
        p.Y = [ ymin, ymax ]
    elseif p.y !== nothing
        p.X = [ xmin, xmax ]
        p.Y = [ p.y, p.y ]
    end
end


function configure!(c::Chart, legend::Legend)
    legend.handle_length = 1.9*legend.fontsize
    legend.row_sep = 0.3*legend.fontsize
    legend.col_sep = 1.5*legend.fontsize
    legend.inner_pad = 1.5*legend.row_sep
    legend.outer_pad = legend.inner_pad

    plots = [ p for p in c.dataseries if p.label != ""]
    
    nlabels = length(plots)
    label_width = maximum( getsize(plot.label, legend.fontsize)[1] for plot in plots )
    label_heigh = maximum( getsize(plot.label, legend.fontsize)[2] for plot in plots )

    handle_length = legend.handle_length
    row_sep       = legend.row_sep
    inner_pad     = legend.inner_pad
    col_sep       = legend.col_sep
    ncols         = legend.ncols

    nrows = ceil(Int, nlabels/legend.ncols)

    col_witdhs = zeros(ncols)
    for k in 1:length(plots)
        j = k%ncols==0 ? ncols : k%ncols # column
        item_width = handle_length + 2*inner_pad + label_width
        col_witdhs[j] = max(col_witdhs[j], item_width)
    end

    legend.height = nrows*label_heigh + (nrows-1)*row_sep + 2*inner_pad
    legend.width = sum(col_witdhs) + (ncols-1)*col_sep + 2*inner_pad

    if legend.location in (:outertopright, :outerright, :outerbottomright)
        c.rightpad += legend.width + c.outerpad
    elseif legend.location in ( :outertopleft, :outerleft, :outerbottomleft)
        c.leftpad += legend.width + c.outerpad
    elseif legend.location == :outertop
        c.toppad += legend.height + 2*c.outerpad
    elseif legend.location == :outerbottom
        c.bottompad += legend.height + 2*c.outerpad
    end

end


function draw!(c::Chart, cc::CairoContext, canvas::Canvas)
    # draw grid
    set_source_rgb(cc, 0.9, 0.9, 0.9) # gray
    set_line_width(cc, 0.2)

    xmin, xmax = c.xaxis.limits
    for x in c.xaxis.ticks
        min(xmax, xmin) <= x<= max(xmax,xmin) || continue
        x1 = canvas.box[1] + c.xaxis.width/(xmax-xmin)*(x-xmin)
        y1 = canvas.box[2]
        y2 = canvas.box[4]
        move_to(cc, x1, y1); line_to(cc, x1, y2); stroke(cc)
    end

    ymin, ymax = c.yaxis.limits
    for y in c.yaxis.ticks
        min(ymax, ymin) <= y<= max(ymax,ymin) || continue
        y1 = canvas.box[2] + c.yaxis.height/(ymax-ymin)*(ymax-y)
        x1 = canvas.box[1]
        x2 = canvas.box[3]
        move_to(cc, x1, y1); line_to(cc, x2, y1); stroke(cc)
    end

    # draw border
    set_source_rgb(cc, 0.0, 0.0, 0.0)
    set_line_width(cc, 0.5)
    x, y = canvas.box[1:2]
    w, h = canvas.box[3:4] - canvas.box[1:2]
    rectangle(cc, x, y, w, h)
    stroke(cc)
end


function draw!(chart::Chart, cc::CairoContext, p::LineSeries)

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
        x1, y1 = data2user(chart.canvas, X[1], Y[1])

        if p.ls==:solid
            move_to(cc, x1, y1)
            for i in 2:n
                x, y = data2user(chart.canvas, X[i], Y[i])
                line_to(cc, x, y);
            end
            stroke(cc)
        else # dashed
            len = sum(p.dash)
            offset = 0.0
            set_dash(cc, p.dash, offset)
            move_to(cc, x1, y1)
            for i in 2:n
                x, y = data2user(chart.canvas, X[i], Y[i])
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
        x, y = data2user(chart.canvas, x, y)
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
        x, y = data2user(chart.canvas, x, y)
        x1, y1 = data2user(chart.canvas, X[i], Y[i])
        x2, y2 = data2user(chart.canvas, X[i+1], Y[i+1])
        α = -atand(y2-y1, x2-x1) # tilt

        # pads
        pad = chart.args.fontsize*0.3
        
        dx = pad*abs(sind(α))
        dy = pad*abs(cosd(α))

        # Default location "top"
        if p.tagloc==:top
            va = "bottom"
            if 0<α<= 90 || -180 <α<= -90 
                ha = "right" 
                dx, dy = -dx, -dy
            else
                ha = "left"
                dy = -dy
            end
        else
            va = "top"
            if 0<α<=90 || -180<α<=-90
                ha = "left"
            else
                ha = "right"
                dx = -dx
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


function draw!(c::Chart, cc::CairoContext, legend::Legend)

    plots = [ p for p in c.dataseries if p.label != ""]
    
    set_font_size(cc, legend.fontsize)
    font = get_font(legend.font)
    select_font_face(cc, font, Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL )

    handle_length = legend.handle_length
    row_sep       = legend.row_sep
    inner_pad     = legend.inner_pad
    outer_pad     = legend.outer_pad
    col_sep       = legend.col_sep
    ncols         = legend.ncols
    
    # update the width    
    col_witdhs = zeros(ncols)
    for (k, plot) in enumerate(plots)
        j = k%ncols==0 ? ncols : k%ncols # column
        label_width = getsize(cc, plot.label, legend.fontsize)[1]
        item_width = handle_length + 2*inner_pad + label_width
        col_witdhs[j] = max(col_witdhs[j], item_width)
    end

    # legend.height = nrows*label_heigh + (nrows-1)*row_sep + 2*inner_pad
    legend.width = sum(col_witdhs) + (ncols-1)*col_sep + 2*inner_pad
    
    # set legend location
    if legend.location in (:topright, :right, :bottomright)
        x1 = c.canvas.box[3] - outer_pad - legend.width
    elseif legend.location in (:top, :bottom, :outertop, :outerbottom)
        x1 = 0.5*(c.canvas.box[1]+c.canvas.box[3]) - legend.width/2
    elseif legend.location in (:topleft, :left, :bottomleft)
        x1 = c.canvas.box[1] + outer_pad
    elseif legend.location in (:outertopleft, :outerleft, :outerbottomleft)
        x1 = c.outerpad
    elseif legend.location in (:outertopright, :outerright, :outerbottomright)
        x1 = c.width - legend.width - c.outerpad
    end
    
    if legend.location in (:topleft, :top, :topright)
        y1 = c.canvas.box[2] + outer_pad
    elseif legend.location in (:left, :right, :outerleft, :outerright)
        y1 = 0.5*(c.canvas.box[2]+c.canvas.box[4]) - legend.height/2
    elseif legend.location in (:bottomleft, :bottom, :bottomright)
        y1 = c.canvas.box[4] - outer_pad - legend.height
    elseif legend.location == :outertop
        y1 = c.outerpad
    elseif legend.location == :outerbottom
        y1 = c.height - legend.height - c.outerpad
    elseif legend.location in (:outertopleft, :outertopright)
        y1 = c.canvas.box[2]
    elseif legend.location in (:outerbottomleft, :outerbottomright)
        y1 = c.canvas.box[4] - legend.height
    end
    x2, y2 = x1 + legend.width, y1 + legend.height

    set_matrix(cc, CairoMatrix([1, 0, 0, 1, 0, 0]...)) 

    # draw rounded rectangle
    r = 0.02*min(c.canvas.width, c.canvas.height)
    move_to(cc, x1, y1+r)
    line_to(cc, x1, y2-r)
    curve_to(cc, x1, y2, x1, y2, x1+r, y2)
    line_to(cc, x2-r, y2)
    curve_to(cc, x2, y2, x2, y2, x2, y2-r)
    line_to(cc, x2, y1+r)
    curve_to(cc, x2, y1, x2, y1, x2-r, y1)
    line_to(cc, x1+r, y1)
    curve_to(cc, x1, y1, x1, y1, x1, y1+r)
    close_path(cc)
    set_source_rgb(cc, 1, 1, 1) # white
    fill_preserve(cc)
    set_source_rgb(cc, 0, 0, 0) # black
    set_line_width(cc, 0.4)
    stroke(cc)

    # draw labels
    label_heigh = maximum( getsize(plot.label, legend.fontsize)[2] for plot in plots )

    for (k, plot) in enumerate(plots)
        i = ceil(Int, k/ncols)  # line
        j = k%ncols==0 ? ncols : k%ncols # column
        x2 = x1 + inner_pad + sum(col_witdhs[1:j-1]) + (j-1)*col_sep

        y2 = y1 + inner_pad + label_heigh/2 + (i-1)*(label_heigh+row_sep)

        set_source_rgb(cc, plot.linecolor...)
        if plot.ls!=:none
            move_to(cc, x2, y2)
            rel_line_to(cc, handle_length, 0)
            set_line_width(cc, plot.lw)
            plot.ls!=:solid && set_dash(cc, plot.dash)
            stroke(cc)
            set_dash(cc, Float64[])
        end

        # draw marker
        x = x2 + handle_length/2
        draw_marker(cc, x, y2, plot.marker, plot.markersize, plot.markercolor, plot.mscolor)

        # draw label
        x = x2 + handle_length + 2*inner_pad
        y = y2 + 0.25*legend.fontsize

        set_source_rgb(cc, 0, 0, 0)
        draw_text(cc, x, y, plot.label, halign="left", valign="bottom", angle=0)
    end

end


function draw!(c::Chart, cc::CairoContext)
    # draw canvas grid
    draw!(c, cc, c.canvas)

    # draw axes
    x = c.leftpad+c.yaxis.width
    y = c.toppad+c.yaxis.height
    move_to(cc, x, y)
    draw!(cc, c.xaxis)
    
    x = c.leftpad
    y = c.toppad
    move_to(cc, x, y)
    draw!(cc, c.yaxis)

    # draw plots
    x, y = c.canvas.box[1:2]
    w, h = c.canvas.box[3:4] - c.canvas.box[1:2]
    rectangle(cc, x, y, w, h)
    Cairo.clip(cc)

    # draw dataseries
    sorted = sort(c.dataseries, by=x->x.order)
    for p in sorted
        draw!(c, cc, p)
    end
    reset_clip(cc)

    # draw annotations
    for a in c.annotations
        draw!(c, cc, a)
    end

    # draw legend
    if c.legend !== nothing
        draw!(c, cc, c.legend)
    end
end
