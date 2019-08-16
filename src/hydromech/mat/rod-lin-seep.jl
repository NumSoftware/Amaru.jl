# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export RodLinSeep

mutable struct RodLinSeepIpState<:IpState
    env::ModelEnv
    V::Float64       # fluid velocity
    uw::Float64      # pore pressure
    function RodLinSeepIpState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.V = 0.0
        this.uw = 0.0
        return this
    end
end

mutable struct RodLinSeep<:Material
    k ::Float64     # specific permeability
    γw::Float64     # specific weight of the fluid
    β ::Float64     # compressibility of fluid
    A ::Float64     # section area

    function RodLinSeep(prms::Dict{Symbol,Float64})
        return  RodLinSeep(;prms...)
    end

    function RodLinSeep(;k=NaN, gammaw=NaN, beta=0.0, A=NaN)
        k<=0.0 && error("Invalid value for k: $k")
        gammaw<=0.0 && error("Invalid value for gammaw: $gammaw")
        beta<0.0 && error("Invalid value for beta: $beta")
        A<=0.0 && error("Invalid value for A: $A")
        return new(k, gammaw, beta, A)
    end
end

matching_elem_type(::RodLinSeep) = SeepRod
#matching_elem_type_if_embedded(::RodLinSeep) = SeepEmbRod

# Create a new instance of Ip data
new_ip_state(mat::RodLinSeep, env::ModelEnv) = RodLinSeepIpState(env)

function update_state!(mat::RodLinSeep, ipd::RodLinSeepIpState, Δuw::Float64, G::Float64)
    k = mat.k
    ipd.V  -= k*G 
    ipd.uw += Δuw
    return ipd.V
end

function ip_state_vals(mat::RodLinSeep, ipd::RodLinSeepIpState)
    return OrderedDict(
      :va => ipd.V,
      :uwa => ipd.uw,
      :Fa => ipd.V*mat.A,
      :A  => mat.A )
end


