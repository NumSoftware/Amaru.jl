# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct Canvas<:ChartComponent
    width::Float64
    height::Float64
    limits::Vector{Float64}
    box::Vector{Float64}
    function Canvas()
        return new(0.0, 0.0)
    end
end


function configure!(c::Chart, canvas::Canvas)
    xmin, xmax = c.xaxis.limits
    ymin, ymax = c.yaxis.limits
    # canvas.limits = [ xmin, ymin, xmax, ymax ]

    if c.args.aspectratio==:equal
        # compute extra limits
        width = c.width - c.yaxis.width - c.outerpad - c.rightpad
        height = c.height - c.xaxis.height - c.toppad - c.outerpad
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
    
    canvas.width = c.width - c.yaxis.width - c.outerpad - c.rightpad
    canvas.height = c.height - c.xaxis.height - c.toppad - c.outerpad
    canvas.box = [ c.outerpad + c.yaxis.width, c.toppad, c.width-c.rightpad, c.height - c.xaxis.height-c.outerpad ]
    canvas.limits = [ xmin, ymin, xmax, ymax ]
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


