# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticLSJoint, ElasticRSJoint

mutable struct ElasticLSJointState<:IpState
    env::ModelEnv
    σ ::Array{Float64,1}
    u ::Array{Float64,1}
    function ElasticLSJointState(env::ModelEnv)
        this = new(env)
        this.σ = zeros(env.ndim)
        this.u = zeros(env.ndim)
        return this
    end
end

ElasticLSJoint_params = [
    FunInfo(:ElasticLSJoint, "Elastic material for a rod-solid interface."),
    KwArgInfo(:ks, "Shear stiffness", cond=:(ks>=0)),
    KwArgInfo(:kn, "Normal stiffness", cond=:(kn>0)),
]
@doc docstring(ElasticLSJoint_params) ElasticLSJoint(; kwargs...)

mutable struct ElasticLSJoint<:Material
    ks::Float64
    kn::Float64

    function ElasticLSJoint(; kwargs...)
        args = checkargs(kwargs, ElasticLSJoint_params)
        this = new(args.ks, args.kn)
        return this
    end
end

const ElasticRSJoint = ElasticLSJoint


# Type of corresponding state structure
compat_state_type(::Type{ElasticLSJoint}, ::Type{MechRSJoint}, env::ModelEnv) = ElasticLSJointState


function calcD(mat::ElasticLSJoint, state::ElasticLSJointState)
    ks = mat.ks
    kn = mat.kn
    if state.env.ndim==2
        return [  ks  0.0
                 0.0   kn ]
    else
        return  [  ks  0.0  0.0
                  0.0   kn  0.0
                  0.0  0.0   kn ]
    end
end


function update_state!(mat::ElasticLSJoint, state::ElasticLSJointState, Δu)
    D = calcD(mat, state)
    Δσ = D*Δu

    state.u .+= Δu
    state.σ .+= Δσ
    return Δσ, success()
end


function ip_state_vals(mat::ElasticLSJoint, state::ElasticLSJointState)
    return OrderedDict(
      :ur   => state.u[1] ,
      :tau  => state.σ[1] )
end
