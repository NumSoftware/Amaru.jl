# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


export Compressible

"""
    $(TYPEDEF)

An elastoplastic constitutive model is used to simulate materials that primarily exhibit
elastic behavior but exhibit nonlinear behavior when the isotropic compression p exceeds 
a predefined yield limit J1y.
"""
mutable struct Compressible<:Material
    E   ::Float64  # initial Young modulus
    ν   ::Float64
    J1y ::Float64
    H   ::Float64

    function Compressible(; args...)
        args_info = arg_rules(Compressible)
        args = checkargs(args, args_info)
        this = new(args.E, args.nu, args.J1y, args.H)
        return this
    end
end

arg_rules(::Type{Compressible}) = 
[
    @arginfo E E>=0.0 "Young modulus"
    @arginfo nu 0.0<=nu<0.5 "Poisson ratio"
    @arginfo J1y J1y<0.0 "Yield first invariant"
    @arginfo H=0.0 H>=0 "Plastic modulus"
]


mutable struct CompressibleState<:IpState
    env ::ModelEnv
    σ   ::Vec6  # current stress
    ε   ::Vec6  # current strain
    εvp ::Float64  # current plastic volumetric strain
    Δλ  ::Float64  # plastic multiplier

    function CompressibleState(env::ModelEnv)
        this      = new(env)
        this.σ    = zeros(Vec6)
        this.ε    = zeros(Vec6)
        this.εvp  = 0.0
        this.Δλ   = 0.0

        return this
    end
end


compat_state_type(::Type{Compressible}, ::Type{MechSolid}, env::ModelEnv) = CompressibleState


@inline function yield_func(mat::Compressible, state::CompressibleState, σ::Vec6,  εvp::Float64)
    return -calcJ1(σ) + mat.J1y + mat.H*εvp
end


function calcD(mat::Compressible, state::CompressibleState)
    De  = calcDe(mat.E, mat.ν, state.env.ana.stressmodel)

    state.Δλ==0.0 && return De

    dfdσ = -I2
    dgdσ = state.σ
    H    = mat.H
    J1   = calcJ1(state.σ)
    
    return De - De*dgdσ*dfdσ'*De/(dfdσ'*De*dgdσ - H*J1)
end


function update_state!(mat::Compressible, state::CompressibleState, Δε::Array{Float64,1})
    σini = state.σ

    De  = calcDe(mat.E, mat.ν, state.env.ana.stressmodel)
    σtr = state.σ + De*Δε
    ftr = yield_func(mat, state, σtr, state.εvp)
    
    if ftr < 1.e-8 # elastic
        state.Δλ = 0.0
        state.σ  = σtr
    else # plastic
        # println()
        # @show ftr
        E, ν, H, J1y = mat.E, mat.ν, mat.H, mat.J1y
        K    = E/(3*(1-2*ν))
        G    = E/(2*(1+ν))
        J1tr = tr(σtr)

        Δλ = -ftr/(3*K*J1y + 3*K*H*state.εvp + H*J1tr )

        J1    = J1tr/(1+3*K*Δλ)
        s     = 1/(1+2*G*Δλ)*dev(σtr)
        state.σ  = s + J1/3*I2

        state.Δλ   = Δλ
        state.εvp += state.Δλ*J1

        # @show state.Δλ
        # @show yield_func(mat, state, state.σ, state.εvp)
        # @assert state.Δλ>0
        # error()
    end


    state.ε += Δε
    Δσ       = state.σ - σini
    return Δσ, success()
end


function ip_state_vals(mat::Compressible, state::CompressibleState)
    dict = stress_strain_dict(state.σ, state.ε, state.env.ana.stressmodel)
    dict[:evp] = state.εvp
    return dict
end
