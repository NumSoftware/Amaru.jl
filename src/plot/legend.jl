# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct Legend<:ChartComponent
    location::Symbol
    font::String
    fontsize::Float64
    handle_length::Float64
    label_sep::Float64
    inner_pad::Float64
    outer_pad::Float64
    width::Float64
    height::Float64

    function Legend(; location=:topright, font="NewComputerModern", fontsize=7)
        return new(location, font, fontsize)
    end
end

function configure!(c::Chart, legend::Legend)
    legend.handle_length = 0.06*c.figsize[1]
    legend.label_sep = 0.3*legend.fontsize
    legend.inner_pad = 1.5*legend.label_sep
    legend.outer_pad = 0.015*c.figsize[1]

    plots = [ p for p in c.dataseries if p.label != ""]
    
    nlabels = length(plots)
    label_width = maximum( getsize(plot.label, legend.fontsize)[1] for plot in plots )
    label_heigh = maximum( getsize(plot.label, legend.fontsize)[2] for plot in plots )

    handle_length = legend.handle_length
    label_sep = legend.label_sep
    inner_pad = legend.inner_pad

    legend.height = nlabels*label_heigh + (nlabels-1)*label_sep + 2*inner_pad 
    legend.width = handle_length + 2*inner_pad + label_width + 2*inner_pad
end

function draw!(c::Chart, cc::CairoContext, legend::Legend)

    set_font_size(cc, legend.fontsize)

    font = get_font(legend.font)
    select_font_face(cc, font, Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL )

    handle_length = legend.handle_length
    label_sep = legend.label_sep
    inner_pad = legend.inner_pad
    outer_pad = legend.outer_pad

    # set legend location
    if legend.location in (:topright, :right, :bottomright)
        x1 = c.canvas.box[3] - outer_pad - legend.width
    elseif legend.location in (:top, :bottom)
        x1 = 0.5*(c.canvas.box[1]+c.canvas.box[3]) - legend.width/2
    else
        x1 = c.canvas.box[1] + outer_pad
    end
    
    if legend.location in (:topleft, :top, :topright)
        y1 = c.canvas.box[2] + outer_pad
    elseif legend.location in (:left, :right)
        y1 = 0.5*(c.canvas.box[2]+c.canvas.box[4]) - legend.height/2
    else
        y1 = c.canvas.box[4] - outer_pad - legend.height
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
    plots = [ p for p in c.dataseries if p.label != ""]
    label_heigh = maximum( getsize(plot.label, legend.fontsize)[2] for plot in plots )
    for (i,plot) in enumerate(plots)
        # draw line
        x = x1 + inner_pad
        y = y1 + inner_pad + label_heigh/2 + (i-1)*(label_heigh+label_sep) 
        move_to(cc, x, y)
        rel_line_to(cc, handle_length, 0)
        set_source_rgb(cc, plot.lc...)
        set_line_width(cc, plot.lw)
        if plot.ls!=:solid
            set_dash(cc, plot.dash)
        end
        stroke(cc)
        set_dash(cc, Float64[])

        x = x1 + inner_pad + handle_length/2
        draw_marker(cc, x, y, plot.ms, plot.marker)

        x = x1 + inner_pad + handle_length + 2*inner_pad
        y = y - 0.15*legend.fontsize
        set_source_rgb(cc, 0, 0, 0)
        draw_text(cc, x, y, plot.label, halign="left", valign="center", angle=0)
    end

end