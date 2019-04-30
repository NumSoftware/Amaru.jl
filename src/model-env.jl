# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct ModelEnv
    ndim::Int              # Analysis dimension
    modeltype::Symbol      # Analysis type (eg. :general, :plane_stress, :plane_strain, etc.
    thickness::Float64     # Model thickness
    transient::Bool        # Time dependent analysis
    t::Float64             # Time in time dependent analysis
    cstage::Int            # Current stage
    cinc::Int              # Current increment
    cout::Int              # Current output file number
    function ModelEnv()
        this = new()
        this.ndim = 3
        this.modeltype = :general # plane_stress, axisymmetric
        this.thickness = 1.0
        this.transient = false
        this.t = 0.0

        this.cstage = 0 # current stage
        this.cinc   = 0 # current increment
        this.cout   = 0 # output files counter
        return this
    end
end
