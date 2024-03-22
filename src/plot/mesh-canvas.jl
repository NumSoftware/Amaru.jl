# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct MeshCanvas<:ChartComponent
    width::Float64
    height::Float64
    limits::Vector{Float64}
    box::Vector{Float64}
    function MeshCanvas()
        return new(0.0, 0.0)
    end
end

function configure!(chart, canvas::MeshCanvas)
end