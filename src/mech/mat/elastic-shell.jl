# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticShell

mutable struct ElasticShellIpState<:IpState
    "Environment information"
    env::ModelEnv
    "Stress tensor"
    σ::Array{Float64,1} # in the local system
    "Strain tensor"
    ε::Array{Float64,1} # in the local system

    function ElasticShellIpState(env::ModelEnv=ModelEnv())
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

    function ElasticShell(prms::Dict{Symbol,Float64})
        return  ElasticShell(;prms...)
    end

    function ElasticShell(;E=NaN, nu=NaN, thickness=0.1)
        E>0.0 || error("Invalid value for E: $E")
        (0<=nu<0.5) || error("Invalid value for nu: $nu")
        thickness >0.0 || error("Invalid value for thickness: $thickness")
        this = new(E, nu, thickness)
        return this
    end
end

matching_elem_type(::ElasticShell) = MechShell

# Type of corresponding state structure
ip_state_type(::ElasticShell) = ElasticShellIpState


function calcD(mat::ElasticShell, ipd::ElasticShellIpState)
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

function stress_update(mat::ElasticShell, ipd::ElasticShellIpState, dε::Array{Float64,1})
    D = calcD(mat, ipd)
    dσ = D*dε
    ipd.ε += dε
    ipd.σ += dσ
    return dσ, success()
end

function ip_state_vals(mat::ElasticShell, ipd::ElasticShellIpState)
    return stress_strain_dict(ipd.σ, ipd.ε, ipd.env.modeltype)
end