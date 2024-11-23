# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LinearCohesive

mutable struct LinearCohesiveState<:IpState
    ctx::Context
    σ   ::Array{Float64,1}
    w   ::Array{Float64,1}
    h   ::Float64
    function LinearCohesiveState(ctx::Context)
        this = new(ctx)
        this.σ = zeros(3)
        this.w = zeros(3)
        this.h = 0.0
        return this
    end
end

LinearCohesive_params = [
    FunInfo(:LinearCohesive, "Consitutive model for cohesive elements with linear elastic behavior."),
    KwArgInfo(:E, "Young modulus", cond=:(E>0)),
    KwArgInfo(:nu, "Poisson ratio", cond=:(0<=nu<0.5)),
    KwArgInfo(:zeta, "Elastic displacement scale factor", cond=:(zeta>0)),
    ArgOpt((:E,:nu), (:kn,:ks))
]
@doc docstring(LinearCohesive_params) LinearCohesive(; kwargs...)

mutable struct LinearCohesive<:Material
    E::Float64 # Young modulus from bulk material
    ν::Float64 # Poisson ratio from bulk material
    ζ::Float64  # elastic displacement scale factor (formerly α)

    function LinearCohesive(; kwargs...)
        args = checkargs(kwargs, LinearCohesive_params)
        this = new(args.E, args.nu, args.zeta)
        return this
    end

end


# Type of corresponding state structure
compat_state_type(::Type{LinearCohesive}, ::Type{MechJoint}, ctx::Context) = LinearCohesiveState


function calcD(mat::LinearCohesive, state::LinearCohesiveState)
    ndim = state.ctx.ndim
    G    = mat.E/(1.0+mat.ν)/2.0
    kn   = mat.E*mat.ζ/state.h
    ks   = G*mat.ζ/state.h

    if ndim==2
        return [  kn  0.0
                 0.0   ks ]
    else
        return  [  kn  0.0  0.0
                  0.0   ks  0.0
                  0.0  0.0   ks ]
    end
end


function update_state!(mat::LinearCohesive, state::LinearCohesiveState, Δu)
    ndim = state.ctx.ndim
    D  = calcD(mat, state)
    Δσ = D*Δu

    state.w[1:ndim] += Δu
    state.σ[1:ndim] += Δσ
    return Δσ, success()
end


function ip_state_vals(mat::LinearCohesive, state::LinearCohesiveState)
    ndim = state.ctx.ndim
    if ndim == 3
       return Dict(
          :jw1  => state.w[1],
          :js1  => state.σ[1],
          )
    else
        return Dict(
          :jw1  => state.w[1],
          :js1  => state.σ[1],
          )
    end
end


function output_keys(mat::LinearCohesive)
    return Symbol[:jw1, :jw1]
end