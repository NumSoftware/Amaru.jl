# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticSolidLinCond

mutable struct ElasticSolidLinCondIpState<:IpState
    env::ModelEnv
    σ::Array{Float64,1}
    ε::Array{Float64,1}
    V::Array{Float64,1}
    uw::Float64
    function ElasticSolidLinCondIpState(env::ModelEnv=ModelEnv()) 
        this = new(env)
        this.σ = zeros(6)
        this.ε = zeros(6)
        this.V = zeros(env.ndim)
        this.uw = 0.0
        return this
    end
end


mutable struct ElasticSolidLinCond<:Material
    E ::Float64
    nu::Float64
    k ::Float64 # thermal conductivity
    α ::Float64 # thermal expansion coefficient

    function ElasticSolidLinCond(prms::Dict{Symbol,Float64})
        return  ElasticSolidLinCond(;prms...)
    end

    function ElasticSolidLinCond(;E=1.0, nu=0.0, k=NaN, gammaw=NaN, alpha=1.0)
        E<=0.0       && error("Invalid value for E: $E")
        !(0<=nu<0.5) && error("Invalid value for nu: $nu")
        isnan(k)     && error("Missing value for k")
        0<=alpha<=1  || error("Invalid value for alpha: $alpha")
        this = new(E, nu, k, gammaw, alpha)
        return this
    end
end

# Returns the element type that works with this material model
matching_elem_type(::ElasticSolidLinCond) = HMSolid

# Create a new instance of Ip data
new_ip_state(mat::ElasticSolidLinCond, env::ModelEnv) = ElasticSolidLinCondIpState(env)

function calcD(mat::ElasticSolidLinCond, ipd::ElasticSolidLinCondIpState)
    return calcDe(mat.E, mat.nu, ipd.env.modeltype) # function calcDe defined at elastic-solid.jl
end

function calcK(mat::ElasticSolidLinCond, ipd::ElasticSolidLinCondIpState) # Hydraulic conductivity matrix
    if ipd.env.ndim==2
        return mat.k*eye(2)
    else
        return mat.k*eye(3)
    end
end

function stress_update(mat::ElasticSolidLinCond, ipd::ElasticSolidLinCondIpState, Δε::Array{Float64,1}, Δuθ::Float64, G::Array{Float64,1})
    De = calcD(mat, ipd)
    Δσ = De*Δε
    ipd.ε  += Δε
    ipd.σ  += Δσ
    K = calcK(mat, ipd)
    ipd.V   = -K*G
    ipd.uw += Δuw
    return Δσ, ipd.V
end

function ip_state_vals(mat::ElasticSolidLinCond, ipd::ElasticSolidLinCondIpState)
    D = stress_strain_dict(ipd.σ, ipd.ε, ipd.env.ndim)

    D[:vx] = ipd.V[1]
    D[:vy] = ipd.V[2]
    if ipd.env.ndim==3
        D[:vz] = ipd.V[3]
    end

    return D
end
