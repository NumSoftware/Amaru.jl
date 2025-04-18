# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export DruckerPrager

DruckerPrager_params = [
    FunInfo(:DruckerPrager, "Linear-elastic model with Drucker-Prager yield criterion and linear hardening"),
    KwArgInfo(:E, "Young modulus", cond=:(E>0.0)),
    KwArgInfo(:nu, "Poisson ratio", 0.0, cond=:(0.0<=nu<0.5)),
    KwArgInfo(:alpha, "alpha parameter", cond=:(alpha>=0.0)),
    KwArgInfo(:kappa, "kappa parameter", cond=:(kappa>=0.0)),
    KwArgInfo(:H, "Hardening modulus", 0.0, cond=:(H>=0.0)),
    KwArgInfo(:rho, "Density", 0.0, cond=:(rho>=0.0))
]
@doc docstring(DruckerPrager_params) DruckerPrager
mutable struct DruckerPrager<:Material
    E::Float64
    ν::Float64
    α::Float64
    κ::Float64
    H::Float64
    ρ::Float64

    function DruckerPrager(;E=NaN, nu=0.0, alpha=0.0, kappa=0.0, H=0.0, rho=0.0)
        @assert E>0.0
        @assert 0.0<=nu<0.5
        @assert alpha>=0.0
        @assert kappa>0.0
        @assert H>=0.0
        @assert rho>=0.0

        this    = new(E, nu, alpha, kappa, H, rho)
        return this
    end
end


"""
    DruckerPragerState

A type for the state data of a `DruckerPrager` type.

# Fields

$(TYPEDFIELDS)
"""
mutable struct DruckerPragerState<:IpState
    ctx::Context
    σ::Vec6
    ε::Vec6
    εpa::Float64
    Δγ::Float64
    function DruckerPragerState(ctx::Context)
        this = new(ctx)
        this.σ   = zeros(Vec6)
        this.ε   = zeros(Vec6)
        this.εpa = 0.0
        this.Δγ  = 0.0
        this
    end
end


# Type of corresponding state structure
compat_state_type(::Type{DruckerPrager}, ::Type{MechSolid}, ctx::Context) = DruckerPragerState

# Element types that work with this material
# compat_elem_types(::Type{DruckerPrager}) = (MechSolid,)


function nlE(fc::Float64, εc::Float64, ε::Array{Float64,1})
    εv = abs(sum(ε[1:3]))
    return 2*fc*(εc-εv)/εc^2
end


function yield_func(mat::DruckerPrager, state::DruckerPragerState, σ::AbstractArray)
    j1  = tr(σ)
    j2d = J2(σ)
    α,κ = mat.α, mat.κ
    H   = mat.H
    εpa = state.εpa
    return α*j1 + √j2d - κ - H*εpa
end


function calcD(mat::DruckerPrager, state::DruckerPragerState)
    α   = mat.α
    H   = mat.H
    De  = calcDe(mat.E, mat.ν, state.ctx.stressmodel)

    if state.Δγ==0.0
        return De
    end

    j2d = J2(state.σ)
    if j2d != 0.0
        s  = dev(state.σ)
        su = s/norm(s)
        V  = α*I2 + su/√2 # df/dσ
        N  = V
        Nu = N/norm(N)
    else # apex
        Nu = 1.0/√3.0*I2
        V  = Nu
    end

    # return De - inner(De,Nu) ⊗ inner(V,De) / (inner(V,De,Nu) + H)
    return De - De*Nu*V'*De / (V'*De*Nu + H)
end


function update_state!(mat::DruckerPrager, state::DruckerPragerState, Δε::Array{Float64,1})
    σini = state.σ
    De   = calcDe(mat.E, mat.ν, state.ctx.stressmodel)
    σtr  = state.σ + De*Δε
    ftr  = yield_func(mat, state, σtr)

    if ftr < 1.e-8
        # elastic
        state.Δγ = 0.0
        state.σ  = σtr
    else
        # plastic
        K, G  = mat.E/(3.0*(1.0-2.0*mat.ν)), mat.E/(2.0*(1.0+mat.ν))
        α, H  = mat.α, mat.H
        n     = 1.0/√(3.0*α*α+0.5)
        j1tr  = tr(σtr)
        j2dtr = J2(σtr)

        if √j2dtr - state.Δγ*n*G > 0.0 # conventional return # todo fix this
            state.Δγ = ftr/(9*α*α*n*K + n*G + H)
            j1     = j1tr - 9*state.Δγ*α*n*K
            m      = 1.0 - state.Δγ*n*G/√j2dtr
            state.σ  = m*dev(σtr) + j1/3.0*I2
        else # return to apex
            κ      = mat.κ
            state.Δγ = (α*j1tr-κ-H*state.εpa)/(3*√3*α*K + H)
            j1     = j1tr - 3*√3*state.Δγ*K
            state.σ  = j1/3.0*I2
        end

        state.εpa += state.Δγ

    end

    state.ε += Δε
    Δσ     = state.σ - σini
    return Δσ, success()
end


function ip_state_vals(mat::DruckerPrager, state::DruckerPragerState)
    σ, ε  = state.σ, state.ε
    j1    = tr(σ)
    srj2d = √J2(σ)

    D = stress_strain_dict(σ, ε, state.ctx.stressmodel)
    D[:epa]   = state.εpa
    D[:j1]    = j1
    D[:srj2d] = srj2d

    return D
end
