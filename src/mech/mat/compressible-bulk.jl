# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export CompressibleBulk

"""
    $(TYPEDEF)

An elastoplastic constitutive model is used to simulate materials that primarily exhibit
elastic behavior but exhibit nonlinear behavior when the isotropic compression p exceeds 
a predefined yield limit py.
"""
mutable struct CompressibleBulk<:Material
    E   ::Float64  # initial Young modulus
    ν   ::Float64
    py0 ::Float64
    pc  ::Float64
    εvpc::Float64 # ≈ 0.07
    α   ::Float64 # ≈ 0.07
    ρ   ::Float64

    function CompressibleBulk(; E=NaN, nu=0.2, pc=NaN, alpha=2.0, epsy=NaN, epsc=NaN, rho=0.0)
        @check E>0.0  
        @check nu>=0.0 
        @check pc<0.0     # pc ≈ 20*fc ≈ 7*py0; py0 ≈ fc/3
        @check epsy<0.0   # epsy ≈ 0.002
        @check epsc<epsy  # epsy ≈ 0.07
        @check alpha>1.0
        @check rho>=0.0

        K = E/(3*(1-2*nu))
        py0 = K*epsy
        εvpc = epsc-epsy

        this = new(E, nu, py0, pc, εvpc, alpha, rho)
        return this
    end
end


mutable struct CompressibleBulkState<:IpState
    env ::ModelEnv
    σ   ::Array{Float64,1}  # current stress
    ε   ::Array{Float64,1}  # current strain
    εvp ::Float64           # current plastic volumetric strain
    Δλ  ::Float64          # plastic multiplier

    function CompressibleBulkState(env::ModelEnv=ModelEnv())
        this      = new(env)
        this.σ    = zeros(6)
        this.ε    = zeros(6)
        this.εvp  = 0.0
        this.Δλ   = 0.0

        return this
    end
end


# Returns the element type that works with this material model
matching_elem_type(::CompressibleBulk) = MechSolid

# Type of corresponding state structure
ip_state_type(::CompressibleBulk) = CompressibleBulkState


@inline function calc_py(mat::CompressibleBulk, state::CompressibleBulkState, εvp::Float64)
    return mat.py0 + (mat.pc-mat.py0)*exp( 1 - (εvp/mat.εvpc)^-mat.α )
end


@inline function yield_func(mat::CompressibleBulk, state::CompressibleBulkState, σ::Array{Float64,1},  εvp::Float64)
    p = (σ[1]+σ[2]+σ[3])/3
    f = calc_py(mat, state, εvp) - p
    @show f
    return calc_py(mat, state, εvp) - p
end


function calcD(mat::CompressibleBulk, state::CompressibleBulkState)
    α   = mat.α
    De  = calcDe(mat.E, mat.ν, state.env.modeltype)

    state.Δλ==0.0 && return De

    εvp  = state.εvp
    εvpc = mat.εvpc
    α    = mat.α

    dpydεvp = (pc-py0)*exp( 1 - (εvp/εvpc)^-α )*α*(εvp/εvpc)^(-α-1)/εvpc

    return De - De*dfdσ*dfdσ*D/(dot(dfdσ*De, dfdσ) + dpydεvp) # todo
end

function calc_σ_Δεvp_Δλ(mat::CompressibleBulk, state::CompressibleBulkState, σtr::Array{Float64,1})
    Δλ = 0.0
    up = 0.0
    σ  = zeros(ndim)
    σ0 = zeros(ndim)
    K  = mat.E/(3*(1-2*mat.ν))
    
    tol    = 1e-6
    maxits = 50
    for i in 1:maxits
        # quantities at n+1
        σ = σtr - Δλ*K*mI

        dfdσ  = -mI/3
        r     = dfdσ
        Δεvp  = Δλ*r
        εvp   = state.εvp + Δεvp

        f       = yield_func(mat, state, σ, εvp)
        dpydεvp = (pc-py0)*exp( 1 - (εvp/εvpc)^-α )*α*(εvp/εvpc)^(-α-1)/εvpc
        dfdΔλ   = dpydεvp - 3*K
        Δλ      = Δλ - f/dfdΔλ
        
        if Δλ<=0 || isnan(Δλ) || i==maxits
            return 0.0, state.σ, 0.0, failure("TCJoint: failed to find Δλ")
        end

        maximum(abs, σ-σ0) <= tol && break
        σ0 .= σ
    end

    return σ, up, Δλ, success()
end


function stress_update(mat::CompressibleBulk, state::CompressibleBulkState, Δε::Array{Float64,1})
    σini = state.σ

    De  = calcDe(mat.E, mat.ν, state.env.modeltype)
    σtr = state.σ + inner(De, Δε)
    ftr = yield_func(mat, state, σtr, state.εvp)

    if ftr < 1.e-8 # elastic
        state.Δλ = 0.0
        state.σ  = σtr
    else
        # plastic
        σ, Δεvp, Δλ, status = calc_σ_Δεvp_Δλ(mat::AbstractTCJoint, state::AbstractTCJointState, σtr::Array{Float64,1})
        failed(status) && return state.σ, status

        state.σ, state.Δεvp, state.Δλ = σ, Δεvp, Δλ

        state.εvp += -Δλ

    end

    @show Δε

    state.ε += Δε
    Δσ       = state.σ - σini
    return Δσ, success()
end


function ip_state_vals(mat::CompressibleBulk, state::CompressibleBulkState)
    dict = stress_strain_dict(state.σ, state.ε, state.env.modeltype)
    dict[:evp] = state.εvp
    return dict
end
