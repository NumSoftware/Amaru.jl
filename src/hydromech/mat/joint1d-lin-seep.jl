# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export Joint1DLinSeep

mutable struct Joint1DLinSeepIpState<:IpState
    env::ModelEnv
    ndim::Int
    V::Float64       # fluid velocity
    uw::Float64      # pore pressure
    function Joint1DLinSeepIpState(env::ModelEnv=ModelEnv())
        this = new(env)
        ndim = env.ndim
        this.V = 0.0
        this.uw = 0.0
        return this
    end
end

mutable struct Joint1DLinSeep<:Material
    k ::Float64    # specific permeability
    γw::Float64    # specific weight of the fluid
    h ::Float64    # section perimeter

    function Joint1DLinSeep(prms::Dict{Symbol,Float64})
        return  Joint1DLinSeep(;prms...)
    end

    function Joint1DLinSeep(;k=NaN, gammaw=NaN, h=NaN, A=NaN, dm=NaN)
        # A : section area
        # dm: section diameter
        # h : section perimeter
        k>=0.0 || error("Invalid value for k: $k")
        gammaw>=0.0 || error("Invalid value for gammaw: $gammaw")
        (h>0 || A>0 || dm>0) || error("perimeter h, section area A or diameter dm should be provided")

        if isnan(h) 
            if A>0
                h = 2.0*(A*pi)^0.5
            else
                h = pi*dm
            end
        end
        @assert h>0

        this = new(k, gammaw, h)
        return this
    end
end


# Returns the element type that works with this material
matching_elem_type(::Joint1DLinSeep) = SeepJoint1D

# Create a new instance of Ip data
new_ip_state(mat::Joint1DLinSeep, env::ModelEnv) = Joint1DLinSeepIpState(env)

function update_state!(mat::Joint1DLinSeep, ipd::Joint1DLinSeepIpState, Δuw::Float64, G::Float64)
    k = mat.k
    ipd.V  -= k*G 
    ipd.uw += Δuw
    return ipd.V
end

function ip_state_vals(mat::Joint1DLinSeep, ipd::Joint1DLinSeepIpState)
    return OrderedDict(
      :va => ipd.V,
      :uwa => ipd.uw)
end
