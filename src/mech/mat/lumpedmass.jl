# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LumpedMass

mutable struct LumpedMassState<:IpState
    env::ModelEnv
    function LumpedMassState(env::ModelEnv=ModelEnv())
        return new(env)
    end
end

mutable struct LumpedMass<:MatParams
    m::Float64

    function LumpedMass(prms::Dict{Symbol,Float64})
        return  LumpedMass(;prms...)
    end

    function LumpedMass(;m::Real=NaN)
        m>=0.0 || error("Invalid value for m: $m")
        return new(m)
    end
end

matching_elem_type(::LumpedMass) = MechLumpedMassElem

# Type of corresponding state structure
ip_state_type(matparams::LumpedMass) = LumpedMassState


function ip_state_vals(matparams::LumpedMass, state::LumpedMassState)
    return OrderedDict{Symbol, Float64}()
end
