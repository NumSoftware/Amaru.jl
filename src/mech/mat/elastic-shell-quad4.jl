# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticShellQUAD4

mutable struct ElasticShellQuad4State<:IpState
    env::ModelEnv
    σ::Array{Float64,1}
    function ElasticShellQuad4State(env::ModelEnv=ModelEnv())
        σ = zeros(5)
        return new(env, σ)
    end
end

mutable struct ElasticShellQUAD4<:Material
    E::Float64
    ν::Float64
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



# Type of corresponding state structure
ip_state_type(::Type{ElasticShellQUAD4}) = ElasticShellQuad4State


function ip_state_vals(mat::ElasticShellQUAD4, state::ElasticShellQuad4State)
    return OrderedDict{Symbol, Float64}()
end
