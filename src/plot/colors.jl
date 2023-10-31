# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

const _default_colors = [ :C1, :C2, :C3, :C4, :C5, :C6, :C7, :C8, :C9, :C10, :C11, :C12, :C13, :C14, :C15, :C16, ]

const _colors_dict = Dict(
    :C1          => (0.0,0.605,0.978),
    :C2          => (0.888,0.435,0.278),
    :C3          => (0.242,0.643,0.304),
    :C4          => (0.764,0.444,0.824),
    :C5          => (0.675,0.555,0.094),
    :C6          => (0.0,0.66575,0.68099),
    :C7          => (0.930,0.36747,0.575),
    :C8          => (0.776,0.50974,0.146),
    :C9          => (0.0,0.66426,0.55295),
    :C10         => (0.558,0.593,0.117),
    :C11         => (0.0,0.66087,0.79817),
    :C12         => (0.609,0.499,0.911),
    :C13         => (0.380,0.551,0.966),
    :C14         => (0.942,0.375,0.451),
    :C15         => (0.868,0.395,0.713),
    :C16         => (0.423,0.622,0.198),
    :aliceblue   => (0.941, 0.973, 1.0),
    :blue        => (0.0, 0.0, 1.0),
    :black       => (0.0, 0.0, 0.0),
    :brown       => (0.647, 0.165, 0.165),
    :cadetblue   => (0.373, 0.62, 0.627),
    :cyan        => (0.0, 1.0, 1.0),
    :darkblue    => (0.0, 0.0, 0.545),
    :darkgreen   => (0.0, 0.392, 0.0),
    :darkmagenta => (0.545, 0.0, 0.545),
    :darkgray    => (0.663, 0.663, 0.663),
    :darkorange  => (1.0, 0.549, 0.0),
    :darkred     => (0.545, 0.0, 0.0),
    :darkcyan    => (0.0, 0.545, 0.545),
    :gray        => (0.502, 0.502, 0.502),
    :green       => (0.0, 0.502, 0.0),
    :grey        => (0.502, 0.502, 0.502),
    :indianred   => (0.804, 0.361, 0.361),
    :indigo      => (0.294, 0.0, 0.51),
    :lightblue   => (0.678, 0.847, 0.902),
    :lightgreen  => (0.565, 0.933, 0.565),
    :magenta     => (1.0, 0.0, 1.0),
    :olive       => (0.502, 0.502, 0.0),
    :orange      => (1.0, 0.647, 0.0),
    :pink        => (1.0, 0.753, 0.796),
    :purple      => (0.502, 0.0, 0.502),
    :red         => (1.0, 0.0, 0.0),
    :royalblue   => (0.255, 0.412, 0.882),
    :steelblue   => (0.275, 0.51, 0.706),
    :violet      => (0.933, 0.51, 0.933),
    :white       => (1.0, 1.0, 1.0),
    :yellow      => (1.0, 1.0, 0.0),
)

const _colors_list = collect(keys(_colors_dict))


struct Colormap
    stops::Vector{Float64}
    colors::Array{Vec3}

    function Colormap(stops, colors)
        @assert length(stops)==length(colors)
        # @assert stops[1]==0.0
        # @assert stops[end]==1.0
        return new(stops, colors)
    end
end

# Interpolate a color
function (cmap::Colormap)(rval)
    rval<=cmap.stops[1] && return cmap.colors[1]
    rval>=cmap.stops[end] && return cmap.colors[end]
    idx = findfirst(>(rval), cmap.stops)

    t = (rval-cmap.stops[idx-1])/(cmap.stops[idx]-cmap.stops[idx-1])
    return (1-t)*cmap.colors[idx-1] + t*cmap.colors[idx]
end

import Amaru.resize
function Amaru.resize(cmap, min, max; diverging=false)
    rmin = cmap.stops[1]
    rmax = cmap.stops[end]

    (min>=0 || max<=0) && (diverging=false)

    if diverging
        rmid = 0.5*(rmax-rmin)

        stops = []
        for rval in cmap.stops
            if rval<rmid
                s = min - min*(rval-rmin)/(rmid-rmin)
            else
                s = max*(rval-rmid)/(rmax-rmid)
            end
            push!(stops, s)
        end
    else
        stops = [ min + (rval-rmin)/(rmax-rmin)*(max-min) for rval in cmap.stops ]
    end

    return Colormap(stops, cmap.colors)
end


const coolwarm = Colormap(
    [ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 ],
    [
        Vec3(0.23, 0.299, 0.754),
        Vec3(0.348, 0.466, 0.888),
        Vec3(0.484, 0.622, 0.975),
        Vec3(0.619, 0.744, 0.999),
        Vec3(0.754, 0.83, 0.961),
        Vec3(0.867, 0.864, 0.863),
        Vec3(0.947, 0.795, 0.717),
        Vec3(0.968, 0.674, 0.557),
        Vec3(0.932, 0.519, 0.406),
        Vec3(0.839, 0.322, 0.265),
        Vec3(0.706, 0.016, 0.15),
    ]
)