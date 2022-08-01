# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticSolidThermo

mutable struct ElasticSolidThermoState<:IpState
    env::ModelEnv
    σ::Array{Float64,1} # stress
    ε::Array{Float64,1} # strain
    QQ::Array{Float64,1} # heat flux
    D::Array{Float64,1}
    ut::Float64
    function ElasticSolidThermoState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ = zeros(6)
        this.ε = zeros(6)
        this.QQ = zeros(env.ndim)
        this.D = zeros(env.ndim) 
        this.ut = 0.0
        return this
    end
end


mutable struct ElasticSolidThermo<:Material
    E ::Float64 # Young's Modulus kN/m2
    nu::Float64 # Poisson coefficient
    k ::Float64 # thermal conductivity  w/m/k
    ρ ::Float64 # density Ton/m3
    cv::Float64 # Specific heat J/Ton/k
    α ::Float64 # thermal expansion coefficient  1/K or 1/°C


    function ElasticSolidThermo(prms::Dict{Symbol,Float64})
        return  ElasticSolidThermo(;prms...)
    end

    function ElasticSolidThermo(;E=NaN, nu=NaN, k=NaN, rho=NaN, cv=NaN, alpha=1.0)
        E>0.0       || error("Invalid value for E: $E")
        (0<=nu<0.5) || error("Invalid value for nu: $nu")
        isnan(k)     && error("Missing value for k")
        cv<=0.0       && error("Invalid value for cv: $cv")
        0<=alpha<=1  || error("Invalid value for alpha: $alpha")

        this = new(E, nu, k, rho, cv, alpha)

        return this
    end
end

# Returns the element type that works with this material model
matching_elem_type(::ElasticSolidThermo) = TMSolid

# Type of corresponding state structure
ip_state_type(mat::ElasticSolidThermo) = ElasticSolidThermoState


#function set_state(ipd::ElasticSolidLinSeepState; sig=zeros(0), eps=zeros(0))
#    sq2 = √2.0
#    mdl = [1, 1, 1, sq2, sq2, sq2]
#    if length(sig)==6
#        ipd.σ .= sig.*mdl
#    else
#        if length(sig)!=0; error("ElasticSolidLinSeep: Wrong size for stress array: $sig") end
#    end
#    if length(eps)==6
#        ipd.ε .= eps.*mdl
#    else
#        if length(eps)!=0; error("ElasticSolidLinSeep: Wrong size for strain array: $eps") end
#    end
#end



function calcD(mat::ElasticSolidThermo, ipd::ElasticSolidThermoState)
    return calcDe(mat.E, mat.nu, ipd.env.modeltype) # function calcDe defined at elastic-solid.jl
end

function calcK(mat::ElasticSolidThermo, ipd::ElasticSolidThermoState) # Thermal conductivity matrix
    if ipd.env.ndim==2
        return mat.k*eye(2)
    else
        return mat.k*eye(3)
    end
end

function stress_update(mat::ElasticSolidThermo, ipd::ElasticSolidThermoState, Δε::Array{Float64,1}, Δut::Float64, G::Array{Float64,1}, Δt::Float64)
    De = calcD(mat, ipd)
    Δσ = De*Δε
    ipd.ε  += Δε
    ipd.σ  += Δσ
    K = calcK(mat, ipd)
    ipd.QQ = -K*G
    ipd.D  += ipd.QQ*Δt
    ipd.ut += Δut
    return Δσ, ipd.QQ
end

function ip_state_vals(mat::ElasticSolidThermo, ipd::ElasticSolidThermoState)
    D = stress_strain_dict(ipd.σ, ipd.ε, ipd.env.modeltype)

    #D[:qx] = ipd.QQ[1] # VERIFICAR NECESSIDADE
    #D[:qy] = ipd.QQ[2] # VERIFICAR NECESSIDADE
    #if ipd.env.ndim==3 # VERIFICAR NECESSIDADE
        #D[:qz] = ipd.QQ[3] # VERIFICAR NECESSIDADE
    #end # VERIFICAR NECESSIDADE

    return D
end
