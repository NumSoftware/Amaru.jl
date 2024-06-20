# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LinearElasticFluid

LinearElasticFluid_params = [
    FunInfo(:LinearElasticFluid, "Linear-elastic-flow material model"),
    KwArgInfo(:K, "Bulk modulus", cond=:(K>0.0)),
]
@doc docstring(LinearElasticFluid_params) LinearElasticFluid

mutable struct LinearElasticFluid<:Material
    K ::Float64

    function LinearElasticFluid(; kwargs...)
        args = checkargs(kwargs, LinearElasticFluid_params)
        return new(args.K)
    end
end


mutable struct ElasticFluidState<:IpState
    env::ModelEnv
    p::Float64   # hydrostatic pressure
    εv::Float64  # volumetric strain

    function ElasticFluidState(env::ModelEnv)
        this = new(env)
        this.p = 0.0
        this.εv = 0.0
        return this
    end
end


compat_state_type(::Type{LinearElasticFluid}, ::Type{MechFluid}, env::ModelEnv)  = ElasticFluidState


function update_state!(mat::LinearElasticFluid, state::ElasticFluidState, dεv::Float64)
    dp = mat.K*dεv
    state.εv += dεv
    state.p += dp
    return dp, success()
end


function ip_state_vals(mat::LinearElasticFluid, state::ElasticFluidState)
    return OrderedDict(
      :p  => state.p,
      :ev => state.εv
    )
end
