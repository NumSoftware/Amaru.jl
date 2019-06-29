# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LinSeep

mutable struct LinSeepIpState<:IpState
    env::ModelEnv
    V::Array{Float64,1} # fluid velocity
    uw::Float64         # pore pressure
    function LinSeepIpState(env::ModelEnv=ModelEnv()) 
        this = new(env)
        this.uw = 0.0
        this.V  = zeros(env.ndim) 
        return this
    end
end


mutable struct LinSeep<:Material
    k ::Float64 # specific permeability
    γw::Float64 # specific weight of the fluid
    S ::Float64 # Storativity coefficient

    function LinSeep(prms::Dict{Symbol,Float64})
        return  LinSeep(;prms...)
    end

    function LinSeep(;k=NaN, gammaw=NaN, S=0.0)

        k>0         || error("Invalid value for k: $k")
        gammaw>0    || error("Invalid value for gammaw: $gammaw")
        S>=0.0      || error("Invalid value for S: $S")

        this    = new(k, gammaw, S)
        return this
    end
end

# Returns the element type that works with this material model
matching_elem_type(::LinSeep) = SeepSolid

# Create a new instance of Ip data
new_ip_state(mat::LinSeep, env::ModelEnv) = LinSeepIpState(env)

function calcK(mat::LinSeep, ipd::LinSeepIpState) # Hydraulic conductivity matrix
    if ipd.env.ndim==2
        return mat.k*eye(2)
    else
        return mat.k*eye(3)
    end
end

function update_state!(mat::LinSeep, ipd::LinSeepIpState, Δuw::Float64, G::Array{Float64,1})
    K = calcK(mat, ipd)
    ipd.V   = -K*G
    ipd.uw += Δuw
    return ipd.V
end


function ip_state_vals(mat::LinSeep, ipd::LinSeepIpState)
    D = Dict{Symbol, Float64}()
    D[:vx] = ipd.V[1]
    D[:vy] = ipd.V[2]
    if ipd.env.ndim==3
        D[:vz] = ipd.V[3]
    end

    return D
end
