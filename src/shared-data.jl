# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct SharedAnalysisData
    ndim::Int              # Analysis dimension
    model_type::Symbol     # Analysis type (eg. :general, :plane_stress, :plane_strain, etc.
    thickness::Float64     # Model thickness
    t::Float64             # Time in time dependent analysis
    function SharedAnalysisData()
        this = new()
        ndim = 3
        model_type = :general
        t = 0.0
        return this
    end
end

