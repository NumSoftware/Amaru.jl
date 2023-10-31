# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct GeometryCanvas<:ChartComponent
    width::Float64
    height::Float64
    limits::Vector{Float64}
    box::Vector{Float64}
    function GeometryCanvas()
        return new(0.0, 0.0)
    end
end

function configure!(chart::Chart, canvas::GeometryCanvas)
end