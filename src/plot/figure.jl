abstract type Figure 
    # width::Float64
    # height::Float64
end

abstract type AbstractChart<:Figure end
abstract type FigureComponent end
abstract type DataSeries end

const _available_formats = [
    ".pdf",
    ".png",
    ".svg",
]


function save(figure::Figure, filename::String, copypath::String="")
    width, height = figure.width, figure.height

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
    configure!(figure)

    if fmt==".png"
        set_source_rgb(cc, 1.0, 1.0, 1.0) # RGB values for white
        paint(cc)
    end

    draw!(figure, cc)
    
    if fmt==".png"
        write_to_png(surf, filename)
    else
        finish(surf)
    end

    copypath!="" && cp(filename, joinpath(dirname(copypath), basename(filename)), force=true)

    return nothing
end
