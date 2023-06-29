# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticShell

mutable struct ElasticShellState<:IpState
    "Environment information"
    env::ModelEnv
    "Stress tensor"
    σ::Array{Float64,1} # in the local system
    "Strain tensor"
    ε::Array{Float64,1} # in the local system

    function ElasticShellState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ = zeros(6)
        this.ε = zeros(6)
        return this
    end
end

mutable struct ElasticShell<:Material
    E::Float64
    nu::Float64
    th::Float64
    ρ::Float64

    function ElasticShell(prms::Dict{Symbol,Float64})
        return  ElasticShell(;prms...)
    end

    function ElasticShell(;E=NaN, nu=NaN, thickness=NaN, rho=0.0)
        E>0.0 || error("Invalid value for E: $E")
        (0<=nu<0.5) || error("Invalid value for nu: $nu")
        thickness >0.0 || error("Invalid value for thickness: $thickness")
        rho<0.0      && error("Invalid value for rho: $rho")
        this = new(E, nu, thickness, rho)
        return this
    end
end

matching_elem_type(::ElasticShell, shape::CellShape, ndim::Int) = MechShell

# Type of corresponding state structure
ip_state_type(::ElasticShell) = ElasticShellState


function calcD(mat::ElasticShell, state::ElasticShellState)
    E = mat.E
    ν = mat.nu
    c = E/(1.0-ν^2)
    g = E/(1+ν)
    return [
        c    c*ν   0.0  0.0    0.0    0.0
        c*ν  c     0.0  0.0    0.0    0.0
        0.0  0.0   0.0  0.0    0.0    0.0
        0.0  0.0   0.0  5/6*g  0.0    0.0
        0.0  0.0   0.0  0.0    5/6*g  0.0
        0.0  0.0   0.0  0.0    0.0    g ]
    # ezz = -ν/E*(sxx+syy)
end

function stress_update(mat::ElasticShell, state::ElasticShellState, dε::Array{Float64,1})
    D = calcD(mat, state)
    dσ = D*dε
    state.ε += dε
    state.σ += dσ
    return dσ, success()
end

function ip_state_vals(mat::ElasticShell, state::ElasticShellState)
    return stress_strain_dict(state.σ, state.ε, state.env.modeltype)
end