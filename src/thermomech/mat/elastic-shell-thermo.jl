# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticShellThermo

mutable struct ElasticShellThermoState<:IpState
    env::ModelEnv
    σ::Vec6 # stress
    ε::Vec6 # strain
    QQ::Array{Float64,1} # heat flux
    D::Array{Float64,1}
    ut::Float64
    function ElasticShellThermoState(env::ModelEnv)
        this = new(env)
        this.σ = zeros(Vec6)
        this.ε = zeros(Vec6)
        this.QQ = zeros(env.ndim)
        this.D = zeros(env.ndim)
        this.ut = 0.0
        return this
    end
end


ElasticShellThermo_params = [
    FunInfo(:ElasticShellThermo, "A material model for elastic shells with thermal effects"),
    KwArgInfo(:E, "Young modulus", cond=:(E>0.0)),
    KwArgInfo(:nu, "Poisson ratio", 0.0, cond=:(0.0<=nu<0.5)),
    KwArgInfo(:k, "Conductivity", cond=:(k>=0)),
    KwArgInfo(:cv, "Specific heat", 0.0, cond=:(cv>=0)),
    KwArgInfo(:alpha, "Thermal expansion coefficient", cond=:(alpha>=0))
]
@doc docstring(ElasticShellThermo_params) ElasticShellThermo

mutable struct ElasticShellThermo<:Material
    E ::Float64 # Young's Modulus kN/m2
    ν::Float64 # Poisson coefficient
    k ::Float64 # thermal conductivity  w/m/k
    cv::Float64
    α ::Float64 # thermal expansion coefficient  1/K or 1/°C

    function ElasticShellThermo(; args...)
        args = checkargs(args, ElasticShellThermo_params)
        
        return new(args.E, args.nu, args.k, args.cv, args.alpha)
    end
end


# Type of corresponding state structure
compat_state_type(::Type{ElasticShellThermo}, ::Type{TMShell}, env::ModelEnv) = ElasticShellThermoState


function calc_cv(mat::ElasticShellThermo, ut::Float64) # Specific heat
    return mat.cv
end

function calc_α(mat::ElasticShellThermo, ut::Float64)
    return mat.α
end

function calcD(mat::ElasticShellThermo, state::ElasticShellThermoState)
    return calcDe(mat.E, mat.ν, "plane-stress")
end


function calcK(mat::ElasticShellThermo, state::ElasticShellThermoState) # Thermal conductivity matrix
    if state.env.ndim==2
        return mat.k*eye(2)
    else
        return mat.k*eye(3)
    end
end


function update_state!(mat::ElasticShellThermo, state::ElasticShellThermoState, Δε::Array{Float64,1}, Δut::Float64, G::Array{Float64,1}, Δt::Float64)
    De = calcDe(mat.E, mat.ν, "plane-stress")
    Δσ = De*Δε
    state.ε  += Δε
    state.σ  += Δσ

    K = calcK(mat, state)
    state.QQ = -K*G
    state.D  += state.QQ*Δt
    state.ut += Δut
    return Δσ, state.QQ, success()
end


function ip_state_vals(mat::ElasticShellThermo, state::ElasticShellThermoState)
    D = stress_strain_dict(state.σ, state.ε, state.env.ana.stressmodel)
    return D
end
