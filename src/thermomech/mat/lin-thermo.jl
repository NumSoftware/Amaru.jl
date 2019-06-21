# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LinThermo

mutable struct LinThermoIpState<:IpState
    env::ModelEnv
    ut::Float64
    Q::Array{Float64,1}
    function LinThermoIpState(env::ModelEnv=ModelEnv()) 
        this = new(env)
        this.ut = 0.0
        this.Q  = zeros(env.ndim)
        return this
    end
end


mutable struct LinThermo<:Material
    k ::Float64 # Thermal conductivity
    ρ ::Float64 # density
    cv::Float64 # Specific heat

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

# Create a new instance of Ip data
new_ip_state(mat::LinThermo, env::ModelEnv) = LinThermoIpState(env)

function calcK(mat::LinThermo, ipd::LinThermoIpState) # Hydraulic conductivity matrix
    ndim = ipd.env.ndim
    return mat.k*eye(ndim)
end


function update_state!(mat::LinThermo, ipd::LinThermoIpState, Δut::Float64, G::Array{Float64,1})
    K = calcK(mat, ipd)
    ipd.Q   = -K*G
    ipd.ut += Δut
    return ipd.Q
end


function ip_state_vals(mat::LinThermo, ipd::LinThermoIpState)
    D = Dict{Symbol, Float64}()
    D[:qx] = ipd.Q[1]
    D[:qy] = ipd.Q[2]
    if ipd.env.ndim==3
        D[:qz] = ipd.Q[3]
    end
    return D
end
