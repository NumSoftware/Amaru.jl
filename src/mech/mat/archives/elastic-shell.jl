# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LinearElastic

mutable struct ElasticShellState<:IpState
    "Environment information"
    env::ModelEnv
    "Stress tensor"
    σ::Array{Float64,1} # in the local system
    "Strain tensor"
    ε::Array{Float64,1} # in the local system

    function ElasticShellState(env::ModelEnv)
        this = new(env)
        this.σ = zeros(6)
        this.ε = zeros(6)
        return this
    end
end

mutable struct LinearElastic<:Material
    E::Float64
    ν::Float64
    th::Float64
    ρ::Float64

    function LinearElastic(prms::Dict{Symbol,Float64})
        return  LinearElastic(;prms...)
    end

    function LinearElastic(;E=NaN, nu=NaN, thickness=NaN, rho=0.0)
        E>0.0 || error("Invalid value for E: $E")
        (0<=nu<0.5) || error("Invalid value for nu: $nu")
        thickness >0.0 || error("Invalid value for thickness: $thickness")
        rho<0.0      && error("Invalid value for rho: $rho")
        this = new(E, nu, thickness, rho)
        return this
    end
end


# Type of corresponding state structure
compat_state_type(::Type{LinearElastic}) = ElasticShellState

# Element types that work with this material
compat_elem_types(::Type{LinearElastic}) = (MechShell,)


function calcD(mat::LinearElastic, state::ElasticShellState)
    E = mat.E
    ν = mat.ν
    c = E/(1-ν^2)
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


function update_state!(mat::LinearElastic, state::ElasticShellState, dε::Array{Float64,1})
    D = calcD(mat, state)
    dσ = D*dε
    state.ε += dε
    state.σ += dσ
    return dσ, success()
end


function ip_state_vals(mat::LinearElastic, state::ElasticShellState)
    return stress_strain_dict(state.σ, state.ε, state.env.ana.stressmodel)
end