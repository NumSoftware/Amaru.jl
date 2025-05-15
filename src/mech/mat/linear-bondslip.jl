# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LinearLSInterface, ElasticRSJoint, ElasticLSJoint

mutable struct LinearLSInterfaceState<:IpState
    ctx::Context
    σ ::Array{Float64,1}
    u ::Array{Float64,1}
    function LinearLSInterfaceState(ctx::Context)
        this = new(ctx)
        this.σ = zeros(ctx.ndim)
        this.u = zeros(ctx.ndim)
        return this
    end
end

LinearLSInterface_params = [
    FunInfo(:LinearLSInterface, "Elastic material for a rod-solid interface."),
    KwArgInfo(:ks, "Shear stiffness", cond=:(ks>=0)),
    KwArgInfo(:kn, "Normal stiffness", cond=:(kn>0)),
]
@doc docstring(LinearLSInterface_params) LinearLSInterface(; kwargs...)

mutable struct LinearLSInterface<:Material
    ks::Float64
    kn::Float64

    function LinearLSInterface(; kwargs...)
        args = checkargs(kwargs, LinearLSInterface_params)
        this = new(args.ks, args.kn)
        return this
    end
end

const ElasticLSJoint = LinearLSInterface
const ElasticRSJoint = LinearLSInterface


# Type of corresponding state structure
compat_state_type(::Type{LinearLSInterface}, ::Type{MechRSJoint}, ctx::Context) = LinearLSInterfaceState


function calcD(mat::LinearLSInterface, state::LinearLSInterfaceState)
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


function update_state!(mat::LinearLSInterface, state::LinearLSInterfaceState, Δu)
    D = calcD(mat, state)
    Δσ = D*Δu

    state.u .+= Δu
    state.σ .+= Δσ
    return Δσ, success()
end


function ip_state_vals(mat::LinearLSInterface, state::LinearLSInterfaceState)
    return OrderedDict(
      :ur   => state.u[1] ,
      :tau  => state.σ[1] )
end
