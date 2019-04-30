# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LumpedMass

mutable struct LumpedMassIpState<:IpState
    env::ModelEnv
    function LumpedMassIpState(env::ModelEnv=ModelEnv())
        return new(env)
    end
end

mutable struct LumpedMass<:Material
    m::Float64

    function LumpedMass(prms::Dict{Symbol,Float64})
        return  LumpedMass(;prms...)
    end

    function LumpedMass(;m::Real=NaN)
        m>=0.0 || error("Invalid value for m: $m")
        return new(m)
    end
end

matching_elem_type(::LumpedMass) = MechLumpedMass

# Create a new instance of Ip data
new_ip_state(mat::LumpedMass, env::ModelEnv) = LumpedMassIpState(env)

function ip_state_vals(mat::LumpedMass, ipd::LumpedMassIpState)
    return OrderedDict{Symbol, Float64}()
end
