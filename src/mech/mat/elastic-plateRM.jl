# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# Rectangular Reissner Mindlin Plate FEM

export ElasticPlateRM

mutable struct ElasticPlateRMState<:IpState
    env::ModelEnv
    function ElasticPlateRMState(env::ModelEnv=ModelEnv())
        return new(env)
    end
end

mutable struct ElasticPlateRM<:MatParams
    E::Float64
    nu::Float64
    ρ::Float64

    function ElasticPlateRM(prms::Dict{Symbol,Float64})
        return  ElasticPlateRM(;prms...)
    end

    function ElasticPlateRM(;E=NaN, nu=NaN, ρ=0.0)
        E>0.0 || error("Invalid value for E: $E")
        (0<=nu<0.5) || error("Invalid value for nu: $nu")

        this = new(E, nu, ρ)
        return this
    end
end



# Type of corresponding state structure
ip_state_type(matparams::ElasticPlateRM) = ElasticPlateRMState


function ip_state_vals(matparams::ElasticPlateRM, state::ElasticPlateRMState)
    return OrderedDict{Symbol, Float64}()
end
