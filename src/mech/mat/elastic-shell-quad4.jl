# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticShellQUAD4

mutable struct ElasticShellQuad4IpState<:IpState
    env::ModelEnv
    σ::Array{Float64,1}
    function ElasticShellQuad4IpState(env::ModelEnv=ModelEnv())
        σ = zeros(5)
        return new(env, σ)
    end
end

mutable struct ElasticShellQUAD4<:Material
    E::Float64
    nu::Float64
    t::Float64

    function ElasticShellQUAD4(prms::Dict{Symbol,Float64})
        return  ElasticShellQUAD4(;prms...)
    end

    function ElasticShellQUAD4(;E=NaN, nu=NaN, thickness=0.1)
        E>0.0 || error("Invalid value for E: $E")
        (0<=nu<0.5) || error("Invalid value for nu: $nu")
        thickness >0.0 || error("Invalid value for thickness: $thickness")
        this = new(E, nu, thickness)
        return this
    end
end

matching_elem_type(::ElasticShellQUAD4) = ShellQUAD4

# Type of corresponding state structure
ip_state_type(mat::ElasticShellQUAD4) = ElasticShellQuad4IpState

function ip_state_vals(mat::ElasticShellQUAD4, ipd::ElasticShellQuad4IpState)
    return OrderedDict{Symbol, Float64}()
end
