# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


Legend_params = [
    FunInfo(:Legend, "Creates a `Legend` instance."),
    KwArgInfo((:location, :loc), "Location of the legend", :topright, values=_legend_positions ),
    KwArgInfo(:font, "Name of the font", "NewComputerModern", type=AbstractString),
    KwArgInfo(:fontsize, "Size of the font in dpi", 7.0, cond=:(fontsize>0)),
    KwArgInfo(:ncols, "Number of columns in the legend", 1, cond=:(ncols>0)),
]
@doc docstring(Legend_params) Legend

mutable struct Legend<:ChartComponent
    location::Symbol
    font::String
    fontsize::Float64
    ncols::Int
    handle_length::Float64  # length of the line
    row_sep::Float64      # separation between labels
    col_sep::Float64      # separation between labels
    inner_pad::Float64      # padding inside the legend box
    outer_pad::Float64      # padding outside the legend box
    width::Float64          # width of the legend box
    height::Float64         # height of the legend box

    function Legend(; kwargs...)
        args = checkargs(kwargs, Legend_params)
        this = new(args.location, args.font, args.fontsize, args.ncols)
        return this
    end
end


function addlegend!(c::Chart, legend::Legend) 
    c.legend = legend
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
    for (k, plot) in enumerate(plots)
        i = ceil(Int, k/ncols)  # line
        j = k%ncols==0 ? ncols : k%ncols # column
        item_width = handle_length + 2*inner_pad + label_width
        col_witdhs[j] = max(col_witdhs[j], item_width)
    end

    legend.height = nrows*label_heigh + (nrows-1)*row_sep + 2*inner_pad
    legend.width = sum(col_witdhs) + (ncols-1)*col_sep + 2*inner_pad
    # item_width = handle_length + 2*inner_pad + label_width
    # legend.width = item_width*ncols + col_sep*(ncols-1) + 2*inner_pad

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
    nrows         = ceil(Int, length(plots)/legend.ncols)
    
    # update the width    

    col_witdhs = zeros(ncols)
    for (k, plot) in enumerate(plots)
        i = ceil(Int, k/ncols)  # line
        j = k%ncols==0 ? ncols : k%ncols # column
        label_width = getsize(cc, plot.label, legend.fontsize)[1]
        item_width = handle_length + 2*inner_pad + label_width
        col_witdhs[j] = max(col_witdhs[j], item_width)
    end

    # legend.height = nrows*label_heigh + (nrows-1)*row_sep + 2*inner_pad
    legend.width = sum(col_witdhs) + (ncols-1)*col_sep + 2*inner_pad
    
    
    # label_width  = maximum( getsize(cc, plot.label, legend.fontsize)[1] for plot in plots )
    # item_width   = handle_length + 2*inner_pad + label_width
    # legend.width = item_width*ncols + col_sep*(ncols-1) + 2*inner_pad

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
    @show col_witdhs

    for (k, plot) in enumerate(plots)
        i = ceil(Int, k/ncols)  # line
        j = k%ncols==0 ? ncols : k%ncols # column
        # item_width = handle_length + 2*inner_pad + label_width
        # item_width = col_witdhs[j]
        x2 = x1 + inner_pad + sum(col_witdhs[1:j-1]) + (j-1)*col_sep
        # x2 = x1 + inner_pad + (j-1)*(item_width + col_sep)

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