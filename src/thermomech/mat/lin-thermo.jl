# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LinThermo

mutable struct LinThermoIpState<:IpState
    env::ModelEnv
    ut::Float64
    QQ::Array{Float64,1}
    D::Array{Float64,1} #
    function LinThermoIpState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.ut = 0.0
        this.QQ  = zeros(env.ndim)
        this.D  = zeros(env.ndim)
        return this
    end
end



mutable struct LinThermo<:Material
    k ::Float64 # Thermal conductivity w/m/k
    ρ ::Float64 # density Ton/m3
    cv::Float64 # Specific heat J/Ton/k
    #E ::Float64 # Young's modulus kPa
    #nu::Float64 # Poisson ratio
    #α ::Float64 #  coefficient of thermal expansion 1/K or 1/°C

    function LinThermo(prms::Dict{Symbol,Float64})
        return  LinThermo(;prms...)
    end


    function LinThermo(;k=NaN, rho=NaN, cv=NaN)
        k  >= 0.0 || error("Invalid value for k")
        rho>= 0.0 || error("Invalid value for rho")
        cv >= 0.0 || error("Invalid value for cv")
        this = new(k,rho,cv)
        return this
    end
end

# Returns the element type that works with this material model
matching_elem_type(::LinThermo) = ThermoSolid

# Type of corresponding state structure
ip_state_type(mat::LinThermo) = LinThermoIpState

function calcK(mat::LinThermo, ipd::LinThermoIpState) # Thermal conductivity matrix
    if ipd.env.ndim==2
        return mat.k*eye(2)
    else
        return mat.k*eye(3)
    end
end

function update_state!(mat::LinThermo, ipd::LinThermoIpState, Δut::Float64, G::Array{Float64,1}, Δt::Float64)
    K = calcK(mat, ipd)
    ipd.QQ   = -K*G
    ipd.D  += ipd.QQ*Δt
    ipd.ut += Δut
    return ipd.QQ
end


function ip_state_vals(mat::LinThermo, ipd::LinThermoIpState)
    D = Dict{Symbol, Float64}()
    D[:qx] = ipd.QQ[1]
    D[:qy] = ipd.QQ[2]
    if ipd.env.ndim==3
        D[:qz] = ipd.QQ[3]
    end
    return D
end
