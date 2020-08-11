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


"""
Return a dictionary with conventional stress and stress values
from stress and strain tensors defined in Mandel notation.
"""
@inline function stress_strain_dict(σ::Tensor2, ε::Tensor2, env::ModelEnv)
    if env.ndim==2;
        if env.modeltype=="axisymmetric"
            return OrderedDict{Symbol,Float64}(
              :srr => σ[1],
              :syy => σ[2],
              :stt => σ[3],
              :sry => σ[6]/SR2,
              :err => ε[1],
              :eyy => ε[2],
              :ett => ε[3],
              :ery => ε[6]/SR2,
              )
        else
            return OrderedDict{Symbol,Float64}(
              :sxx => σ[1],
              :syy => σ[2],
              :szz => σ[3],
              :sxy => σ[6]/SR2,
              :exx => ε[1],
              :eyy => ε[2],
              :ezz => ε[3],
              :exy => ε[6]/SR2,
              )
        end
    else
        return OrderedDict{Symbol,Float64}(
          :sxx => σ[1],
          :syy => σ[2],
          :szz => σ[3],
          :syz => σ[4]/SR2,
          :sxz => σ[5]/SR2,
          :sxy => σ[6]/SR2,
          :exx => ε[1],
          :eyy => ε[2],
          :ezz => ε[3],
          :eyz => ε[4]/SR2,
          :exz => ε[5]/SR2,
          :exy => ε[6]/SR2,
          )
    end
end
