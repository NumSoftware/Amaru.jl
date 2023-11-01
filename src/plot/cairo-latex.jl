# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


function cairo_set_ft_font(cc, font)
    font_face = ccall( (:cairo_ft_font_face_create_for_ft_face, Cairo.libcairo),
        Ptr{Cvoid}, (FreeTypeAbstraction.FT_Face, Cint),
        font, 0 )
    
    ccall((:cairo_set_font_face, Cairo.libcairo), Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), cc.ptr, font_face)

    return font_face
end


function cairo_font_face_destroy(font_face)
    ccall(
        (:cairo_font_face_destroy, Cairo.libcairo),
        Cvoid, (Ptr{Cvoid},),
        font_face
    )
end
 

function getsize(str::LaTeXString, fontsize::Float64)
    texelems = generate_tex_elements(str)
    width = maximum([elem[2][1] for elem in texelems], init=0.0) + 0.5
    maxh = maximum([elem[2][2] for elem in texelems], init=0.0) + 0.75
    minh = minimum([elem[2][2] for elem in texelems], init=0.0)
    height = maxh - minh
    ffactor = 96/72 # todo: check this factor
    ffactor = 1
    return width*fontsize*ffactor, height*fontsize*ffactor
end

function getsize(str::AbstractString, fontsize::Float64)
    str = LaTeXString(str)
    return getsize(str, fontsize)
end


# Draw latex string with Cairo
function textext(cc::CairoContext, x, y, str::LaTeXString; halign="center", valign="center", angle=0)
    str = latexstring(str)
    texelems = generate_tex_elements(str)
    
    fmatrix = get_font_matrix(cc)
    fsize = norm([fmatrix.xx, fmatrix.xy])

    width, height = getsize(str, fsize)
    # width = getwidth(cc, texelems)
    # height = getheight(cc, texelems)

    rw = halign=="center" ? 0.5 : halign=="right" ? 1.0 : 0.0
    rh = valign=="center" ? 0.5 : valign=="top" ? 1.0 : 0.0

    xref = x + rw*width
    yref = y - rh*height

    θ = angle*pi/180
    T = [ cos(-θ) -sin(-θ); sin(-θ) cos(-θ)]

    ffactor = 96/72  # converding points to dpi
    ffactor = 1

    # for elem in texelems
        
    #     if elem[1] isa TeXChar && elem[3]<1
    #         coord = elem[2]
    #         elem[2] = elem[2] .+ [0, 0.1]
    #     end
    # end

    for elem in texelems
        x0, y0 = elem[2]

        # todo: check these fix of baseline
        # elem[1] isa TeXChar && (y0 -= 0.43*elem[3])
        # y0 += 0.12

        # elem[1] isa TeXChar && (y0 -= 0.43*elem[3])
        # y0 += 0.12

        # user coordinates
        xi = x + x0*fsize*ffactor
        yi = y - y0*fsize*ffactor

        xi, yi = T*[ xi-xref, yi-yref ] .+ [x, y]
        if elem[1] isa TeXChar
            texchar = elem[1]
            scale = elem[3]
            glyph = texchar.represented_char
            # force_italic = occursin(glyph, "αβγδεζηθικλμνξοπρστυφχψω")
            # fnt = texchar.slanted || force_italic ? "NewComputerModern Italic" : "NewComputerModern Regular"
            # set_font_size(cc, fsize*scale)
            set_font_size(cc, fsize*scale)

            # set_font_face(cc, fnt*" "*string(round(fsize*scale, digits=1)))
            # @show fsize*scale

            move_to(cc, xi, yi)
            # @show fsize
            # @show glyph
            show_text(cc, string(glyph))
            # text(cc, xi, yi, string(glyph), halign="left", valign="bottom", angle=angle)
        elseif elem[1] isa HLine
            hline = elem[1]
            w = hline.width*fsize*ffactor
            move_to(cc, xi, yi); rel_line_to(cc, w*cos(-θ), w*sin(-θ)); stroke(cc)
            set_line_width(cc, hline.thickness*fsize)
        end
    end

    set_font_size(cc, fsize)
    set_font_face(cc, "NewComputerModern $fsize") # for pango text

end