# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticBeam

mutable struct BeamIpState<:IpState
    env::ModelEnv
    #σ::Float64
    #ε::Float64
    function BeamIpState(env::ModelEnv=ModelEnv())
        return new(env)
    end
end

mutable struct ElasticBeam<:Material
    E::Float64
    A::Float64
    I::Float64

    function ElasticBeam(prms::Dict{Symbol,Float64})
        return  ElasticBeam(;prms...)
    end

    function ElasticBeam(;E=NaN, A=NaN, I=NaN, Ix=NaN, Iy=NaN, Iz=NaN)
        @assert E>0.0
        @assert A>0.0
        @assert I>0.0
        this = new(E,A,I)
        return this
    end
end

matching_elem_type(::ElasticBeam) = MechBeam

# Create a new instance of Ip data
new_ip_state(mat::ElasticBeam, env::ModelEnv) = BeamIpState(env)

function ip_state_vals(mat::ElasticBeam, ipd::BeamIpState)
    return OrderedDict{Symbol, Float64}()
end
