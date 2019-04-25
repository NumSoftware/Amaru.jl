# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct ModelEnv
    ndim::Int              # Analysis dimension
    modeltype::Symbol      # Analysis type (eg. :general, :plane_stress, :plane_strain, etc.
    thickness::Float64     # Model thickness
    t::Float64             # Time in time dependent analysis
    nstage::Int            # Current stage
    ninc::Int              # Current increment
    nout::Int              # Current output file number
    function ModelEnv()
        this = new()
        this.ndim = 3
        this.modeltype = :general # plane_stress, axisymmetric
        this.thickness = 1.0
        this.t = 0.0

        this.nstage = 0 # current stage
        this.ninc   = 0 # current increment
        this.nout   = 0 # output files counter
        return this
    end
end
