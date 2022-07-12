# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticCookShell

mutable struct ElasticCookShellIpState<:IpState
    "Environment information"
    env::ModelEnv
    "Stress tensor"
    σ::Array{Float64,1}
    "Strain tensor"
    ε::Array{Float64,1}

    function ElasticCookShellIpState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ = zeros(6)
        this.ε = zeros(6)
        return this
    end
end

mutable struct ElasticCookShell<:Material
    E::Float64
    nu::Float64
    t::Float64

    function ElasticCookShell(prms::Dict{Symbol,Float64})
        return  ElasticCookShell(;prms...)
    end

    function ElasticCookShell(;E=NaN, nu=NaN, thickness=0.1)
        E>0.0 || error("Invalid value for E: $E")
        (0<=nu<0.5) || error("Invalid value for nu: $nu")
        thickness >0.0 || error("Invalid value for thickness: $thickness")
        this = new(E, nu, thickness)
        return this
    end
end

matching_elem_type(::ElasticCookShell) = CookShell

# Type of corresponding state structure
ip_state_type(mat::ElasticCookShell) = ElasticCookShellIpState


function calcD(mat::ElasticCookShell, ipd::ElasticCookShellIpState)
    E = mat.E
    ν = mat.nu
    c = E/(1.0-ν^2)
    G = E/(2*(1+ν))
    g = 2*G
    return [
        c    c*ν   0.0  0.0    0.0    0.0
        c*ν  c     0.0  0.0    0.0    0.0
        0.0  0.0   0.0  0.0    0.0    0.0
        0.0  0.0   0.0  5/6*g  0.0    0.0
        0.0  0.0   0.0  0.0    5/6*g  0.0
        0.0  0.0   0.0  0.0    0.0    g ]
    # ezz = -ν/E*(sxx+syy)
end

function stress_update(mat::ElasticCookShell, ipd::ElasticCookShellIpState, dε::Array{Float64,1})
    D = calcD(mat, ipd)
    dσ = D*dε
    ipd.ε += dε
    ipd.σ += dσ
    return dσ, success()
end

function ip_state_vals(mat::ElasticCookShell, ipd::ElasticCookShellIpState)
    return OrderedDict{Symbol,Float64}(
        :sxx => ipd.σ[1],
        :syy => ipd.σ[2],
        :szz => ipd.σ[3],
        :syz => ipd.σ[4],
        :sxz => ipd.σ[5],
        :sxy => ipd.σ[6],
        :exx => ipd.ε[1],
        :eyy => ipd.ε[2],
        :ezz => ipd.ε[3],
        :eyz => ipd.ε[4],
        :exz => ipd.ε[5],
        :exy => ipd.ε[6],
    )
end