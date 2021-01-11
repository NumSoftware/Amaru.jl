# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct ModelEnv
    ndim::Int              # Analysis dimension
    modeltype::String      # Analysis type (eg. "3d", "plane-stress", "plane-strain", etc.
    thickness::Float64     # Model thickness
    transient::Bool        # Time dependent analysis
    t::Float64             # Time in time dependent analysis
    T::Float64             # Pseudo time for current stage
    T0::Float64            # Ambient absolute temperature in Celsius
    cstage::Int            # Current stage
    cinc::Int              # Current increment
    cout::Int              # Current output file number
    outdir::String         # Output directory
    params::Dict{Symbol,Float64} # Global parameters (e.g. gammaw, T0, etc.)
    function ModelEnv()
        this = new()
        this.ndim = 3
        this.modeltype = "3d"
        this.thickness = 1.0
        this.transient = false
        this.t = 0.0
        this.T = 0.0
        this.T0 = 0

        this.cstage = 0 # current stage
        this.cinc   = 0 # current increment
        this.cout   = 0 # output files counter
        this.outdir = ""

        this.params = Dict()
        return this
    end
end
