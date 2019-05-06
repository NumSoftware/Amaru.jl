# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticBeam

mutable struct BeamIpState<:IpState
    env::ModelEnv
    function BeamIpState(env::ModelEnv=ModelEnv())
        return new(env)
    end
end

mutable struct ElasticBeam<:Material
    E::Float64
    A::Float64
    I::Float64
    ρ::Float64

    function ElasticBeam(prms::Dict{Symbol,Float64})
        return  ElasticBeam(;prms...)
    end

    function ElasticBeam(;E=NaN, A=NaN, I=NaN, Ix=NaN, Iy=NaN, Iz=NaN, rho=0.0)
        E>0.0 || error("Invalid value for E: $E")
        A>0.0 || error("Invalid value for A: $A")
        I>0.0 || error("Invalid value for I: $I")
        this = new(E, A, I, rho)
        return this
    end
end

matching_elem_type(::ElasticBeam) = MechBeam

# Create a new instance of Ip data
new_ip_state(mat::ElasticBeam, env::ModelEnv) = BeamIpState(env)

function ip_state_vals(mat::ElasticBeam, ipd::BeamIpState)
    return OrderedDict{Symbol, Float64}()
end
