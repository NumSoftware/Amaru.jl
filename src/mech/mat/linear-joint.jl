# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticJoint

mutable struct JointState<:IpState
    env::ModelEnv
    σ   ::Array{Float64,1}
    w   ::Array{Float64,1}
    h   ::Float64
    function JointState(env::ModelEnv)
        this = new(env)
        this.σ   = zeros(3)
        this.w = zeros(3)
        this.h = 0.0
        return this
    end
end

ElasticJoint_params = [
    FunInfo(:ElasticJoint, "Consitutive model for joints with linear elastic behavior."),
    KwArgInfo(:kn, "Normal stiffness per area", cond=:(kn>0)),
    KwArgInfo(:ks, "Shear stiffness per area", cond=:(ks>=0)),
]
@doc docstring(ElasticJoint_params) ElasticJoint(; kwargs...)

mutable struct ElasticJoint<:Material
    kn::Float64 # Normal stiffness
    ks::Float64 # Shear stiffness

    function ElasticJoint(; kwargs...)
        args = checkargs(kwargs, ElasticJoint_params)
        this = new(args.kn, args.ks)
        return this
    end
end


# Type of corresponding state structure
compat_state_type(::Type{ElasticJoint}, ::Type{MechJoint}, env::ModelEnv) = JointState


function calcD(mat::ElasticJoint, state::JointState)
    ndim = state.env.ndim
    kn   = mat.kn
    ks   = mat.ks

    if ndim==2
        return [  kn  0.0
                 0.0   ks ]
    else
        return  [  kn  0.0  0.0
                  0.0   ks  0.0
                  0.0  0.0   ks ]
    end
end


function update_state!(mat::ElasticJoint, state::JointState, Δu)
    ndim = state.env.ndim
    D  = calcD(mat, state)
    Δσ = D*Δu

    state.w[1:ndim] += Δu
    state.σ[1:ndim] += Δσ
    return Δσ, success()
end


function ip_state_vals(mat::ElasticJoint, state::JointState)
    ndim = state.env.ndim
    if ndim == 3
       return Dict(
          :jw1  => state.w[1],
          :jw2  => state.w[2],
          :jw3  => state.w[3],
          :js1  => state.σ[1],
          :js2  => state.σ[2],
          :js3  => state.σ[3],
          )
    else
        return Dict(
          :jw1  => state.w[1],
          :jw2  => state.w[2],
          :js1  => state.σ[1],
          :js2  => state.σ[2],
          )
    end
end


function output_keys(mat::ElasticJoint)
    return Symbol[:jw1, :jw1]
end