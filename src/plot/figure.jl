abstract type Figure 
    # width::Float64
    # height::Float64
end

abstract type FigureComponent end
abstract type DataSeries end

const _available_formats = [
    ".pdf",
    ".png",
    ".svg",
]


function save(figure::Figure, files::String...)
    configure!(figure)

    for file in files
        width, height = figure.width, figure.height

        fmt = splitext(file)[end]
        if fmt==".pdf"
            surf = CairoPDFSurface(file, width, height)
        elseif fmt==".svg"
            surf = CairoSVGSurface(file, width, height)
        elseif fmt==".ps"
            surf = CairoPSSurface(file, width, height)
        elseif fmt==".png"
            surf = CairoImageSurface(width, height, Cairo.FORMAT_ARGB32)
        else
            formats = join(_available_formats, ", ", " and ")
            throw(AmaruException("Cannot save image to format $fmt. Available formats are: $formats"))
        end

        ctx = CairoContext(surf)  
        set_antialias(ctx, Cairo.ANTIALIAS_NONE) # ANTIALIAS_DEFAULT, ANTIALIAS_NONE, ANTIALIAS_GRAY, ANTIALIAS_SUBPIXEL

        if fmt==".png"
            set_source_rgb(ctx, 1.0, 1.0, 1.0) # RGB values for white
            paint(ctx)
        end

        draw!(figure, ctx)
        
        if fmt==".png"
            write_to_png(surf, file)
        else
            finish(surf)
        end

        figure.quiet || println("  figure saved to $file")

        # copypath!="" && cp(file, joinpath(dirname(copypath), basename(file)), force=true)
    end

end
