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
        # @check alpha>1.0
        @check rho>=0.0

        K = E/(3*(1-2*nu))
        py0 = K*epsy
        εvpc = epsc-epsy
        @check py0>pc

        this = new(E, nu, py0, pc, εvpc, alpha, rho)
        return this
    end
end


mutable struct CompressibleBulkState<:IpState
    env ::ModelEnv
    σ   ::Vec6  # current stress
    ε   ::Vec6  # current strain
    εvp ::Float64           # current plastic volumetric strain
    Δλ  ::Float64          # plastic multiplier

    function CompressibleBulkState(env::ModelEnv=ModelEnv())
        this      = new(env)
        this.σ    = zeros(Vec6)
        this.ε    = zeros(Vec6)
        this.εvp  = 0.0
        this.Δλ   = 0.0

        return this
    end
end


# Type of corresponding state structure
ip_state_type(::MechSolid, ::CompressibleBulk) = CompressibleBulkState


@inline function calc_py(mat::CompressibleBulk, state::CompressibleBulkState, εvp::Float64)
    py0  = mat.py0
    pc   = mat.pc
    α    = mat.α
    εvpc = mat.εvpc
    return py0 + (pc-py0)*exp( 1 - (εvp/εvpc)^-α )
end


@inline function yield_func(mat::CompressibleBulk, state::CompressibleBulkState, σ::Vec6,  εvp::Float64)
    # p = (σ[1]+σ[2]+σ[3])/3
    # f = calc_py(mat, state, εvp) - p
    return calc_py(mat, state, εvp) - tr(σ)/3
end


function calcD(mat::CompressibleBulk, state::CompressibleBulkState)
    α   = mat.α
    De  = calcDe(mat.E, mat.ν, state.env.ana.stressmodel)

    state.Δλ==0.0 && return De

    εvpc = mat.εvpc
    α    = mat.α
    pc   = mat.pc
    py0  = mat.py0
    εvp  = state.εvp

    dfdσ = -tI/3
    dgdσ = dfdσ
    # dgdσ = state.σ
    
    # dpydεvp = (pc-py0)*exp( 1 - (εvp/εvpc)^-α )*α*(εvp/εvpc)^(-α-1)/εvpc
    py = calc_py(mat, state, εvp)
    dpydεvp = (py-py0)*α/εvpc*(εvp/εvpc)^(-α-1)

    return De - De*dgdσ*dfdσ'*De/(dfdσ'*De*dgdσ + dpydεvp)
end


function calc_σ_Δεvp_Δλ(mat::CompressibleBulk, state::CompressibleBulkState, σtr::Vec6)
    Δλ = 0.0
    up = 0.0
    ndim = state.env.ndim
    σ  = zeros(Vec6)
    σ0 = zeros(Vec6)
    K  = mat.E/(3*(1-2*mat.ν))
    pc = mat.pc
    py0 = mat.py0

    εvpc = mat.εvpc
    α = mat.α
    
    tol    = 1e-6
    maxits = 50
    maxits = 20
    Δλ = 1e-11
    εvp = 0.0
    for i in 1:maxits
        # println()
        # @show i
        
        # quantities at n+1
        σ .= σtr .+ Δλ*K*tI
        
        # @show σtr
        # @show σ

        dgdσ = σ
        r    = -tI/3
        Δεp  = Δλ*r
        Δεvp = tr(Δεp)
        εvp  = state.εvp + Δεvp

        f    = yield_func(mat, state, σ, εvp)
        # @show f
        # @show r
        # @show Δεp

        if εvp==0
            dpydεvp = 0.0
        else    
            py = calc_py(mat, state, εvp)
            dpydεvp = (py-py0)*α/εvpc*(εvp/εvpc)^(-α-1)
            # @show εvp
            # @show εvpc
            # @show (εvp/εvpc)^-α
            # @show exp( 1 - (εvp/εvpc)^-α )
            # @show (pc-py0)*exp( 1 - (εvp/εvpc)^-α )
            # @show (py-py0)
            # @show (εvp/εvpc)^(-α-1)
        end

        # @show pc
        # @show py0
        # @show py

        # dfdΔλ   = -dpydεvp - K
        dfdΔλ   = -dpydεvp - K
        Δλ      = Δλ - f/dfdΔλ

        # @show dfdΔλ
        # @show dpydεvp
        # @show Δλ
        
        if Δλ<=0 || isnan(Δλ) || i==maxits
            return state.σ, 0.0, 0.0, failure("TCJoint: failed to find Δλ")
        end

        maximum(abs, σ-σ0) <= tol && break
        σ0 .= σ
    end

    # error()

    return σ, εvp, Δλ, success()
end


function update_state(mat::CompressibleBulk, state::CompressibleBulkState, Δε::Array{Float64,1})
    σini = state.σ

    De  = calcDe(mat.E, mat.ν, state.env.ana.stressmodel)
    σtr = state.σ + De*Δε
    ftr = yield_func(mat, state, σtr, state.εvp)

    if ftr < 1.e-8 # elastic
        state.Δλ = 0.0
        state.σ  = σtr
    else
        # plastic
        σ, εvp, Δλ, status = calc_σ_Δεvp_Δλ(mat, state, σtr)
        failed(status) && return state.σ, status

        state.σ, state.εvp, state.Δλ = σ, εvp, Δλ

        # state.εvp += Δεvp

    end

    # @show Δε

    state.ε += Δε
    Δσ       = state.σ - σini
    return Δσ, success()
end


function ip_state_vals(mat::CompressibleBulk, state::CompressibleBulkState)
    dict = stress_strain_dict(state.σ, state.ε, state.env.ana.stressmodel)
    dict[:evp] = state.εvp
    return dict
end
