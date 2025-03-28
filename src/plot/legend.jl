# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

const _legend_positions=[
    :right,
    :left,
    :top,
    :bottom,
    :topright,
    :topleft,
    :bottomright,
    :bottomleft,
    :outerright,
    :outerleft,
    :outertop,
    :outertopleft,
    :outertopright,
    :outerbottom,
    :outerbottomright,
    :outerbottomleft
]

Legend_params = [
    FunInfo(:Legend, "Creates a `Legend` instance."),
    KwArgInfo((:location, :loc), "Location of the legend", :topright, values=_legend_positions ),
    KwArgInfo(:font, "Name of the font", "NewComputerModern", type=AbstractString),
    KwArgInfo(:fontsize, "Size of the font in dpi", 7.0, cond=:(fontsize>0)),
    KwArgInfo(:ncols, "Number of columns in the legend", 1, cond=:(ncols>0)),
]
@doc docstring(Legend_params) Legend

mutable struct Legend<:FigureComponent
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

