# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

_axes_widget_locations = [:none, :topright, :topleft, :bottomright, :bottomleft]

mutable struct AxisWidget<:FigureComponent
    location::Symbol
    labels::Vector{AbstractString}
    font::String
    font_size::Float64
    azimut::Float64
    elevation::Float64
    arrow_length::Float64
    width::Float64
    height::Float64

    function AxisWidget(; args...)
        args = checkargs(args, 
            [
                KwArgInfo( :location, "AxisWidget location", :bottomleft, values=_axes_widget_locations ),
                KwArgInfo( :labels, "Axis labels", "", type=Vector{AbstractString} ),
                KwArgInfo( :font, "Font name", "NewComputerModern", type=AbstractString),
                KwArgInfo( :font_size, "Font size", 9.0, cond=:(font_size>0)),
                KwArgInfo( :azimut, "Azimut angle for 3d in degrees", 0 ),
                KwArgInfo( :elevation, "Elevation angle for 3d in degrees", 0 ),
                KwArgInfo( :arrow_length, "Length of axis arrow", 20 ),
            ],
            aliens=false,
        )

        this = new(args.location, args.labels, args.font, args.font_size, args.azimut, args.elevation, args.arrow_length)
        return this
    end
end


function configure!(aw::AxisWidget)
    ndim = length(aw.labels)
    if ndim==2
        head_width = 0.15*aw.arrow_length
        aw.width   = aw.arrow_length + head_width/2
        aw.height  = aw.arrow_length + head_width/2
    end
end


function draw_arrow(ctx::CairoContext, x1, y1, x2, y2; head_length=7)
    # Draw the shaft
    move_to(ctx, x1, y1)
    line_to(ctx, x2, y2)
    stroke(ctx)

    θ = π/8

    N  = Vec3(0, 0, 1)
    v1 = normalize(Vec3(x2-x1, y2-y1, 0))
    v3 = v1*cos(π+θ) + cross(N, v1)*sin(π+θ)
    v4 = v1*cos(π-θ) + cross(N, v1)*sin(π-θ)
    p2 = Vec3(x2, y2, 0)
    p0 = p2 - v1*head_length*0.7
    p3 = p2 + v3*head_length/cos(θ/2)
    p4 = p2 + v4*head_length/cos(θ/2)


    move_to(ctx, x2, y2)
    line_to(ctx, p3[1], p3[2])
    line_to(ctx, p0[1], p0[2])
    line_to(ctx, p4[1], p4[2])

    close_path(ctx)
    fill(ctx)

end


function draw!(ctx::CairoContext, aw::AxisWidget)

    Cairo.save(ctx)
    x0, y0 = get_current_point(ctx)
    
    font = get_font(aw.font)
    select_font_face(ctx, font, Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL )
    set_font_size(ctx, aw.font_size)
    set_line_width(ctx, 0.7)

    set_matrix(ctx, CairoMatrix([1, 0, 0, 1, 0, 0]...))
    head_width = 0.15*aw.arrow_length

    translate(ctx, x0+head_width, y0+aw.arrow_length)

    x_label_padding = 0.1*aw.arrow_length
    y_label_padding = 0.18*aw.arrow_length

    set_source_rgb(ctx, _colors_dict[:indianred]...)
    draw_arrow(ctx, 0, 0, aw.arrow_length, 0)
    draw_text(ctx, aw.arrow_length, -x_label_padding, aw.labels[1], halign="right", valign="bottom", angle=0)

    set_source_rgb(ctx, _colors_dict[:green]...)
    draw_arrow(ctx, 0, 0, 0, -aw.arrow_length)
    draw_text(ctx, y_label_padding, -aw.arrow_length, aw.labels[2], halign="left", valign="top", angle=0)


end