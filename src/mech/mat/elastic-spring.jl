# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticSpring

mutable struct ElasticSpringIpState<:IpState
    env::ModelEnv
    function ElasticSpringIpState(env::ModelEnv=ModelEnv())
        return new(env)
    end
end

mutable struct ElasticSpring<:Material
    kx::Float64
    ky::Float64
    kz::Float64
    cx::Float64
    cy::Float64
    cz::Float64

    function ElasticSpring(prms::Dict{Symbol,Float64})
        return  ElasticSpring(;prms...)
    end

    function ElasticSpring(;kx::Number=NaN, ky::Number=0.0, kz::Number=0.0, cx::Number=0.0, cy::Number=0.0, cz::Number=0.0)
        kx>=0.0 || error("Invalid value for kx: $kx")
        ky>=0.0 || error("Invalid value for ky: $ky")
        kz>=0.0 || error("Invalid value for kz: $kz")
        cx>=0.0 || error("Invalid value for cx: $cx")
        cy>=0.0 || error("Invalid value for cy: $cy")
        cz>=0.0 || error("Invalid value for cz: $cz")
        return new(kx, ky, kz, cx, cy, cz)
    end
end

matching_elem_type(::ElasticSpring) = MechSpring

# Type of corresponding state structure
ip_state_type(mat::ElasticSpring) = ElasticSpringIpState

function ip_state_vals(mat::ElasticSpring, ipd::ElasticSpringIpState)
    return OrderedDict{Symbol, Float64}()
end
