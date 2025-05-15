# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LinearJoint

mutable struct LinearJointState<:IpState
    ctx::Context
    σ   ::Array{Float64,1}
    w   ::Array{Float64,1}
    h   ::Float64
    function LinearJointState(ctx::Context)
        this = new(ctx)
        this.σ   = zeros(3)
        this.w = zeros(3)
        this.h = 0.0
        return this
    end
end

LinearJoint_params = [
    FunInfo(:LinearJoint, "Consitutive model for joints with linear elastic behavior."),
    KwArgInfo(:kn, "Normal stiffness per area", cond=:(kn>0)),
    KwArgInfo(:ks, "Shear stiffness per area", cond=:(ks>=0)),
]
@doc docstring(LinearJoint_params) LinearJoint(; kwargs...)

mutable struct LinearJoint<:Material
    kn::Float64 # Normal stiffness
    ks::Float64 # Shear stiffness

    function LinearJoint(; kwargs...)
        args = checkargs(kwargs, LinearJoint_params)
        this = new(args.kn, args.ks)
        return this
    end
end


# Type of corresponding state structure
compat_state_type(::Type{LinearJoint}, ::Type{MechJoint}, ::Context) = LinearJointState


function calcD(mat::LinearJoint, state::LinearJointState)
    ndim = state.ctx.ndim
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


function update_state!(mat::LinearJoint, state::LinearJointState, Δu)
    ndim = state.ctx.ndim
    D  = calcD(mat, state)
    Δσ = D*Δu

    state.w[1:ndim] += Δu
    state.σ[1:ndim] += Δσ
    return Δσ, success()
end


function ip_state_vals(::LinearJoint, state::LinearJointState)
    ndim = state.ctx.ndim
    if ndim == 3
       return Dict(
          :jw  => state.w[1],
          :jw2  => state.w[2],
          :jw3  => state.w[3],
          :jσn  => state.σ[1],
          :js2  => state.σ[2],
          :js3  => state.σ[3],
          )
    else
        return Dict(
          :jw  => state.w[1],
          :jw2  => state.w[2],
          :jσn  => state.σ[1],
          :js2  => state.σ[2],
          )
    end
end


function output_keys(::LinearJoint)
    return Symbol[:jw, :jw]
end