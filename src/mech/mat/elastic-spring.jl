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

    function ElasticSpring(prms::Dict{Symbol,Float64})
        return  ElasticSpring(;prms...)
    end

    function ElasticSpring(;kx::Number=NaN, ky::Number=0.0, kz::Number=0.0)
        kx>=0.0 || error("Invalid value for kx: $kx")
        ky>=0.0 || error("Invalid value for ky: $ky")
        kz>=0.0 || error("Invalid value for kz: $kz")
        #kx+ky+kz>0.0 || error("Invalid values for kx, ky or kz")
        return new(kx, ky, kz)
    end
end

matching_elem_type(::ElasticSpring) = MechSpring

# Create a new instance of Ip data
new_ip_state(mat::ElasticSpring, env::ModelEnv) = ElasticSpringIpState(env)

function ip_state_vals(mat::ElasticSpring, ipd::ElasticSpringIpState)
    return OrderedDict{Symbol, Float64}()
end
