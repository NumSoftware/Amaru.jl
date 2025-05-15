# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export CebBondSlip, CebBondSlip

mutable struct CebBondSlipState<:IpState
    ctx::Context
    σ  ::Array{Float64,1}
    u  ::Array{Float64,1}
    τy ::Float64      # max stress
    sy ::Float64      # accumulated relative displacement
    elastic::Bool
    function CebBondSlipState(ctx::Context)
        this = new(ctx)
        ndim = ctx.ndim
        this.σ = zeros(ndim)
        this.u = zeros(ndim)
        this.τy = 0.0
        this.sy = 0.0
        this.elastic = false
        return this
    end
end


CebBondSlip_params = [
    FunInfo(:CebBondSlip, "Consitutive model for a rod-solid interface according to CEB."),
    KwArgInfo(:taumax, "Shear strength", cond=:(taumax>0)),
    KwArgInfo(:taures, "Residual shear stress", cond=:(taures>=0)),
    KwArgInfo(:s1, "Characteristic slip 1", cond=:(s1>0)),
    KwArgInfo(:s2, "Characteristic slip 2", cond=:(s2>0)),
    KwArgInfo(:s3, "Characteristic slip 3", cond=:(s3>0)),
    KwArgInfo(:alpha, "Ascending curvature parameter", 0.4, cond=:(0.0<=alpha<=1.0)),
    KwArgInfo(:beta, "Descending curvature parameter", 1.0, cond=:(0.0<=beta<=1.0)),
    KwArgInfo(:kn, "Normal stiffness", cond=:(kn>0)),
    KwArgInfo(:ks, "Shear stiffness", 0.0, cond=:(ks>=0)),
    ArgCond(:(taumax>taures)),
    ArgCond(:(ks>=taumax/s1)),
]
@doc docstring(CebBondSlip_params) CebBondSlip(; kwargs...)


mutable struct CebBondSlip<:Material
    τmax:: Float64
    τres:: Float64
    s1  :: Float64
    s2  :: Float64
    s3  :: Float64
    α   :: Float64
    β   :: Float64
    ks  :: Float64
    kn  :: Float64

    function CebBondSlip(; kwargs...)
        args = checkargs(kwargs, CebBondSlip_params)

        this = new(args.taumax, args.taures, args.s1, args.s2, args.s3, args.alpha, args.beta, args.ks, args.kn)
        return this
    end
end

const CebRSJoint = CebBondSlip
const CebBondSlip = CebBondSlip


compat_state_type(::Type{CebBondSlip}, ::Type{MechBondSlip}, ctx::Context) = CebBondSlipState


# Type of corresponding state structure
compat_state_type(::Type{CebBondSlip}) = CebBondSlipState

# Element types that work with this material
compat_elem_types(::Type{CebBondSlip}) = (MechBondSlip,)

CEBJoint1D = CebBondSlip #! deprecated
export CEBJoint1D


function Tau(mat::CebBondSlip, sy::Float64)
    if sy<mat.s1
        return mat.τmax*(sy/mat.s1)^mat.α
    elseif sy<mat.s2
        return mat.τmax
    elseif sy<mat.s3
        return mat.τmax - (mat.τmax-mat.τres)*((sy-mat.s2)/(mat.s3-mat.s2))^mat.β
    else
        return mat.τres
    end
end


function deriv(mat::CebBondSlip, state::CebBondSlipState, sy::Float64)
    if sy==0.0
        s1_factor = 0.01
        sy = s1_factor*mat.s1   # to avoid undefined derivative
    end

    if sy<=mat.s1
        return mat.α*mat.τmax/mat.s1*(sy/mat.s1)^(mat.α-1)
    elseif sy<mat.s2
        return mat.ks*1e-3
    elseif sy<mat.s3
        return -mat.β*(mat.τmax-mat.τres)/(mat.s3-mat.s2)*((sy-mat.s2)/(mat.s3-mat.s2))^(mat.β-1)
    else
        return mat.ks*1e-3
    end
end


function calcD(mat::CebBondSlip, state::CebBondSlipState)
    ks = mat.ks

    if !state.elastic
        dτydsy = deriv(mat, state, state.sy)
        ks = dτydsy
    end

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


function yield_func(mat::CebBondSlip, state::CebBondSlipState, τ::Float64)
    return abs(τ) - state.τy
end

function update_state!(mat::CebBondSlip, state::CebBondSlipState, Δu::Vect)
    ks = mat.ks
    kn = mat.kn
    Δs = Δu[1]      # relative displacement
    τini = state.σ[1] # initial shear stress
    τtr  = τini + ks*Δs # elastic trial

    ftr  = yield_func(mat, state, τtr)

    if ftr<0.0
        τ = τtr
        state.elastic = true
    else
        # dτydsy = deriv(mat, state, state.sy)
        # @s ks
        # @s dτydsy
        # Δsy     = (abs(τtr)-state.τy)/(ks+dτydsy)
        # Δsy     = (abs(τtr)-state.τy)/abs(ks+dτydsy)
        Δsy     = (abs(τtr)-state.τy)/ks
        state.sy += Δsy
        state.τy  = Tau(mat, state.sy)
        τ  = state.τy*sign(τtr)
        Δτ = τ - τini
        state.elastic = false
    end

    # calculate Δσ
    Δτ = τ - τini
    Δσ = kn*Δu
    Δσ[1] = Δτ

    # update u and σ
    state.u .+= Δu
    state.σ .+= Δσ

    return Δσ, success()
end


function ip_state_vals(mat::CebBondSlip, state::CebBondSlipState)
    return OrderedDict(
      :ur   => state.u[1] ,
      :tau  => state.σ[1] ,
      )
end

