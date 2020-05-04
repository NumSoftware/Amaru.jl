# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticPlate

mutable struct PlateIpState<:IpState
    env::ModelEnv
    function PlateIpState(env::ModelEnv=ModelEnv())
        return new(env)
    end
end

mutable struct ElasticPlate<:Material
    E::Float64
    nu::Float64
    ρ::Float64

    function ElasticPlate(prms::Dict{Symbol,Float64})
        return  ElasticPlate(;prms...)
    end

    function ElasticPlate(;E=NaN, nu=NaN, ρ=0.0)
        E>0.0 || error("Invalid value for E: $E")
        (0<=nu<0.5) || error("Invalid value for nu: $nu")

        this = new(E, nu, ρ)
        return this
    end
end

matching_elem_type(::ElasticPlate) = PlateMZC

# Type of corresponding state structure
ip_state_type(mat::ElasticPlate) = PlateIpState

function ip_state_vals(mat::ElasticPlate, ipd::PlateIpState)
    return OrderedDict{Symbol, Float64}()
end
