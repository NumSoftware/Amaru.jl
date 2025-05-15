# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LinearBondSlip, ElasticRSJoint, ElasticBondSlip

mutable struct LinearBondSlipState<:IpState
    ctx::Context
    σ ::Array{Float64,1}
    u ::Array{Float64,1}
    function LinearBondSlipState(ctx::Context)
        this = new(ctx)
        this.σ = zeros(ctx.ndim)
        this.u = zeros(ctx.ndim)
        return this
    end
end

LinearBondSlip_params = [
    FunInfo(:LinearBondSlip, "Elastic material for a rod-solid interface."),
    KwArgInfo(:ks, "Shear stiffness", cond=:(ks>=0)),
    KwArgInfo(:kn, "Normal stiffness", cond=:(kn>0)),
]
@doc docstring(LinearBondSlip_params) LinearBondSlip(; kwargs...)

mutable struct LinearBondSlip<:Material
    ks::Float64
    kn::Float64

    function LinearBondSlip(; kwargs...)
        args = checkargs(kwargs, LinearBondSlip_params)
        this = new(args.ks, args.kn)
        return this
    end
end

const ElasticBondSlip = LinearBondSlip
const ElasticRSJoint = LinearBondSlip


# Type of corresponding state structure
compat_state_type(::Type{LinearBondSlip}, ::Type{MechBondSlip}, ctx::Context) = LinearBondSlipState


function calcD(mat::LinearBondSlip, state::LinearBondSlipState)
    ks = mat.ks
    kn = mat.kn
    if state.ctx.ndim==2
        return [  ks  0.0
                 0.0   kn ]
    else
        return  [  ks  0.0  0.0
                  0.0   kn  0.0
                  0.0  0.0   kn ]
    end
end


function update_state!(mat::LinearBondSlip, state::LinearBondSlipState, Δu)
    D = calcD(mat, state)
    Δσ = D*Δu

    state.u .+= Δu
    state.σ .+= Δσ
    return Δσ, success()
end


function ip_state_vals(mat::LinearBondSlip, state::LinearBondSlipState)
    return OrderedDict(
      :ur   => state.u[1] ,
      :tau  => state.σ[1] )
end
