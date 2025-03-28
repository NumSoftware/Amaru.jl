
Annotation_params = [
    FunInfo(:Annotation, "Creates an `Annotation` instance."),
    ArgInfo(:text, "Text to be displayed", type=AbstractString),
    ArgInfo(:x, "x-coordinate of the annotation", cond=:(0<=x<=1), type=Real),
    ArgInfo(:y, "y-coordinate of the annotation", cond=:(0<=y<=1), type=Real),
    KwArgInfo(:textalignment, "Alignment of the text", :auto, values=(:auto, :left, :right, :top, :bottom)),
    KwArgInfo(:target, "Coordinates of the target point of the arrow relative to data", [0.0,0.0], length=2),
    KwArgInfo((:lw, :lineweight), "Line weight", 0.4, cond=:(lw>0)),
    KwArgInfo(:font, "Name of the font", "NewComputerModern", type=AbstractString),
    KwArgInfo(:fontsize, "Size of the font in dpi", 6.0, cond=:(fontsize>0)),
    KwArgInfo(:color, "Color of the text", :default),
]


mutable struct Annotation <: FigureComponent
    text::AbstractString
    x::Float64
    y::Float64
    textalignment::Symbol
    target::Vector
    lw::Float64
    font::String
    fontsize::Float64
    color::Symbol
    function Annotation(text, x, y; kwargs...)
        args = checkargs([text, x, y], kwargs, Annotation_params)
        target = collect(float.(args.target))
        this = new(args.text, args.x, args.y, args.textalignment, target, args.lw, args.font, args.fontsize, args.color)
        return this
    end
end


function addannotation!(c::AbstractChart, a::Annotation)
    push!(c.annotations, a)
end


function draw!(c::AbstractChart, cc::CairoContext, a::Annotation)
    
    set_font_size(cc, a.fontsize)
    font = get_font(a.font)
    select_font_face(cc, font, Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL )

    set_matrix(cc, CairoMatrix([1, 0, 0, 1, 0, 0]...))

    # convert from axes to Cairo coordinates
    x = c.canvas.box[1] + a.x*c.canvas.width
    y = c.canvas.box[2] + (1-a.y)*c.canvas.height
    halign = a.textalignment==:right ? "right" : "left"
    valign = a.textalignment==:top ? "top" : "bottom"
    set_source_rgb(cc, 0, 0, 0)
    draw_text(cc, x, y, a.text, halign=halign, valign=valign, angle=0)

    if a.textalignment==:auto
        a.textalignment = :left
    end

    # draw arrow
    if a.target !== nothing

        # compute text size
        w, h = getsize(cc, a.text, a.fontsize)
        text_outerpad = 0.1*min(w, h)

        if halign=="left"
            x += w/2
        else
            x -= w/2
        end

        if valign=="top"
            y += 0.5*h
        else
            y -= 0.5*h
        end

        w += text_outerpad
        h += text_outerpad

        # target coordinates
        xa, ya = data2user(c.canvas, a.target[1], a.target[2])

        # deltas
        dx = xa - x
        dy = ya - y

        # compute lines
        if abs(dx)>abs(dy) 
            if abs(dy)<h/2
                lines = "-|"
                if dx>0
                    x += w/2 # right
                else
                    x -= w/2 # left
                end
            else # two lines
                lines = "|-"
                if dy>0 # top
                    y += h/2
                else # bottom
                    y -= h/2
                end
            end
        else
            if abs(dx)<w/2
                lines = "|-"
                if dy>0
                    y += h/2 # top
                else
                    y -= h/2 # bottom
                end
            else # two lines
                lines = "-|"
                if dx>0 # right
                    x += w/2
                else # left
                    x -= w/2
                end
            end
        end

        set_source_rgb(cc, get_color(a.color, :black)...)
        set_line_join(cc, Cairo.CAIRO_LINE_JOIN_ROUND)
        set_line_width(cc, a.lw)
        
        # update deltas
        dx = xa - x
        dy = ya - y

        # Draw line 1
        move_to(cc, x, y)
        if lines[1]=='|'
            rel_line_to(cc, 0, dy)
            y += dy
        else
            rel_line_to(cc, dx, 0)
            x += dx
        end

        # Draw line 2
        dx += sign(dx)*a.lw
        move_to(cc, x, y)
        if lines[2]=='|'
            rel_line_to(cc, 0, dy)
            y += dy
        else
            rel_line_to(cc, dx, 0)
            x += dx
        end

        stroke(cc)
    end
end