# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export CompressiveSolid

"""
    $(TYPEDEF)
"""
mutable struct CompressiveSolid<:Material
    E   ::Float64  # initial Young modulus
    ν   ::Float64
    py0 ::Float64
    pc  ::Float64
    εvpc::Float64 # ≈ 0.07
    ρ   ::Float64

    function CompressiveSolid(; E=NaN, nu=0.2, fc=NaN, pc=20*fc, alpha=2.0, epsy=NaN, epsc=NaN, rho=0.0)
        @check E>0.0  
        @check nu>=0.0 
        @check fc<0.0     # py ≈ fc/3
        @check pc<0.0     # pc ≈ 20*fc
        @check epsy<0.0
        @check epsc<epsy
        @check alpha>=0.0
        @check rho>=0.0

        K = E/(3*(1-2*nu))
        py0 = K*epsy
        εvpc = epsc-epsy

        this = new(E, nu, py0, pc, εvpc, rho)
        return this
    end
end


mutable struct CompressiveSolidState<:IpState
    env ::ModelEnv
    σ   ::Array{Float64,1}  # current stress
    ε   ::Array{Float64,1}  # current strain
    εvp ::Float64           # current plastic volumetric strain
    Δλ  ::Float64          # plastic multiplier

    function CompressiveSolidState(env::ModelEnv=ModelEnv())
        this      = new(env)
        this.σ    = zeros(6)
        this.ε    = zeros(6)
        this.εvp  = 0.0
        this.Δλ   = 0.0

        return this
    end
end


# Returns the element type that works with this material model
matching_elem_type(::CompressiveSolid) = MechSolid

# Type of corresponding state structure
ip_state_type(::CompressiveSolid) = CompressiveSolidState


@inline function calc_py(mat::CompressiveSolid, state::CompressiveSolidState, εvp::Float64)
    return mat.py0 + (mat.pc-map.py0)*exp( 1 - (εvp/mat.εvpc)^-mat.α )
end

function yield_func(mat::CompressiveSolid, state::CompressiveSolidState, σ::Array{Float64,1},  εvp::Float64)
    p = trace(σ)/3
    return calc_py(mat, state, εvp) - p
end

@inline function yield_derivs(mat::AbstractTCJoint, state::AbstractTCJointState, σ::Array{Float64,1}, σmax::Float64)
    return [ -1/3, -1/3, -1/3, 0.0, 0.0, 0.0 ]
end

function calcD(mat::DruckerPrager, state::DruckerPragerState)
    α   = mat.α
    De  = calcDe(mat.E, mat.ν, state.env.modeltype)

    state.Δλ==0.0 && return De

    εvp  = state.εvp
    εvpc = mat.εvpc
    α    = mat.α

    dpydεvp = (pc-py0)*exp( 1 - (εvp/εvpc)^-α )*α*(εvp/εvpc)^(-α-1)/εvpc

    return De - De*dfdσ*dfdσ*D/(dot(dfdσ*De, dfdσ) + dpydεvp) # todo
end



function stress_update(mat::DruckerPrager, state::DruckerPragerState, Δε::Array{Float64,1})
    σini = state.σ
    De   = calcDe(mat.E, mat.ν, state.env.modeltype)
    σtr  = state.σ + inner(De, Δε)
    ftr  = yield_func(mat, state, σtr, state.εvp)

    if ftr < 1.e-8 # elastic
        state.Δλ = 0.0
        state.σ  = σtr
    else
        # plastic
        K, G  = mat.E/(3.0*(1.0-2.0*mat.ν)), mat.E/(2.0*(1.0+mat.ν))
        α, H  = mat.α, mat.H
        n     = 1.0/√(3.0*α*α+0.5)
        j1tr  = J1(σtr)
        j2dtr = J2D(σtr)

        if √j2dtr - state.Δλ*n*G > 0.0 # conventional return
            state.Δλ = ftr/(9*α*α*n*K + n*G + H)
            j1     = j1tr - 9*state.Δλ*α*n*K
            m      = 1.0 - state.Δλ*n*G/√j2dtr
            state.σ  = m*dev(σtr) + j1/3.0*tI
        else # return to apex
            κ      = mat.κ
            state.Δλ = (α*j1tr-κ-H*state.εpa)/(3*√3*α*K + H)
            j1     = j1tr - 3*√3*state.Δλ*K
            state.σ  = j1/3.0*tI
        end

        state.εpa += state.Δλ

    end

    state.ε += Δε
    Δσ     = state.σ - σini
    return Δσ, success()
end

function ip_state_vals(mat::CompressiveSolid, state::CompressiveSolidState)
    dict = stress_strain_dict(state.σ, state.ε, state.env.modeltype)
    dict[:evp] = state.εvp
    return dict
end
