# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

function get_font(str)
    font = FreeTypeAbstraction.findfont(str)
    default_fonts = [
        "NewComputerModern Regular",
        "Times New Roman",
        "Cambria",
        "Palatino Linotype",
        "Georgia",
    ]

    if font===nothing
        for fnt in default_fonts
            font = FreeTypeAbstraction.findfont(fnt)
            font!==nothing && break
        end
    end

    family = replace(font.family_name, "Math" => "")
    return family
end


function getsize(str::LaTeXString, fontsize::Float64)
    texelems = generate_tex_elements(str)
    width = maximum([elem[2][1] for elem in texelems], init=0.0) + 0.7
    maxh = maximum([elem[2][2] for elem in texelems], init=0.0) + 0.7
    minh = minimum([elem[2][2] for elem in texelems], init=0.0) - 0.2
    height = maxh - minh
    return width*fontsize, height*fontsize
end


function getsize(str::AbstractString, fontsize::Float64)
    str = LaTeXString(str)
    return getsize(str, fontsize)
end


function draw_text(cc::CairoContext, x, y, str::AbstractString; halign="center", valign="center", angle=0)
    Cairo.save(cc)

    te = text_extents(cc, str)
    width = te[3]
    height = te[4]
    rw = halign=="center" ? 0.5 : halign=="right" ? 1.0 : 0.0
    rh = valign=="center" ? 0.5 : valign=="top" ? 1.0 : 0.0

    translate(cc, x, y)
    Cairo.rotate(cc, -angle*pi/180)
    translate(cc, -rw*width, rh*height)
    move_to(cc, 0, 0)
    show_text(cc, str)

    Cairo.restore(cc)
end


# Draw latex string with Cairo
function draw_text(cc::CairoContext, x, y, str::LaTeXString; halign="center", valign="center", angle=0)
    Cairo.save(cc)

    texelems = generate_tex_elements(str)
    
    # get the current font
    font_face = ccall((:cairo_get_font_face, Cairo.libcairo), Ptr{Cvoid}, (Ptr{Cvoid},), cc.ptr)
    font_family = ccall((:cairo_toy_font_face_get_family, Cairo.libcairo), Cstring, (Ptr{Cvoid},), font_face)
    font = unsafe_string(font_family)

    # get font size
    fmatrix = get_font_matrix(cc)
    fsize = norm([fmatrix.xx, fmatrix.xy])

    # get size
    width = fsize*(maximum([elem[2][1] for elem in texelems], init=0.0) + 0.7)
    maxh = maximum([elem[2][2] for elem in texelems], init=0.0) + 0.7
    minh = minimum([elem[2][2] for elem in texelems], init=0.0) - 0.15
    height = fsize*(maxh - minh)

    # width, height = getsize(str, fsize)
    rw = halign=="center" ? 0.5 : halign=="right" ? 1.0 : 0.0
    rh = valign=="center" ? 0.5 : valign=="top" ? 1.0 : 0.0

    translate(cc, x, y)
    Cairo.rotate(cc, -angle*pi/180)
    Cairo.translate(cc, -rw*width, rh*height + minh*fsize)
    move_to(cc, 0, 0)


    for elem in texelems
        x0, y0 = elem[2]

        if elem[1] isa TeXChar
            texchar = elem[1]
            scale = elem[3]
            glyph = texchar.represented_char

            if texchar.slanted || occursin(glyph, "αβγδεζηθικλμνξοπρστυφχψω")
                select_font_face(cc, font, Cairo.FONT_SLANT_ITALIC, Cairo.FONT_WEIGHT_NORMAL )
            else
                select_font_face(cc, font, Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL )
            end
            set_font_size(cc, fsize*scale)

            move_to(cc, x0*fsize, -y0*fsize)
            show_text(cc, string(glyph))

        elseif elem[1] isa HLine
            hline = elem[1]
            w = hline.width*fsize
            # move_to(cc, xi, yi); rel_line_to(cc, w*cos(-θ), w*sin(-θ)); stroke(cc)
            set_line_width(cc, hline.thickness*fsize)
            move_to(cc, x0*fsize, -y0*fsize); rel_line_to(cc, w, 0); stroke(cc)
        end
    end

    set_font_size(cc, fsize)
    set_font_face(cc, "NewComputerModern $fsize") # for pango text
    Cairo.restore(cc)

end