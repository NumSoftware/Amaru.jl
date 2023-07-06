# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# mutable struct StageBits
#     T      ::Float64  # Pseudo time for current stage
#     ΔT     ::Float64  # Pseudo time current increment
#     residue::Float64  # Ambient absolute temperature in Celsius
#     stage  ::Int      # Current stage
#     inc    ::Int      # Current increment
#     out    ::Int      # Current output file number
#     inlinearrange::Bool # Flag to chek if analysis is still in linear range
# end

mutable struct ModelEnv
    ndim     ::Int       # Analysis dimension
    transient::Bool      # Time dependent analysis
    t        ::Float64   # Time in time dependent analysis
    outdir   ::String    # Output directory
    anaprops ::AnalysisProps

    # stagebits::StageBits
    T      ::Float64  # Pseudo time for current stage
    ΔT     ::Float64  # Pseudo time current increment
    residue::Float64  # Ambient absolute temperature in Celsius
    stage  ::Int      # Current stage
    inc    ::Int      # Current increment
    out    ::Int      # Current output file number

    function ModelEnv()
        this           = new()
        
        # Analysis related
        this.ndim      = 3
        this.transient = false
        this.t         = 0.0
        this.outdir    = ""
        
        # Stage related:
        # this.stagebits = StageBits(0.0, 0.0, 0.0, 0, 0, 0, true)
        this.T       = 0.0
        this.ΔT      = 0.0
        this.residue = 0.0
        this.stage   = 0
        this.inc     = 0
        this.out     = 0

        return this
    end
end
