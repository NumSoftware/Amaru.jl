# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticShellDegenerate

mutable struct ElasticShellDegenerateIpState<:IpState
    env::ModelEnv
    σ::Array{Float64,1}
    function ElasticShellDegenerateIpState(env::ModelEnv=ModelEnv())
        σ = zeros(5)
        return new(env, σ)
    end
end

mutable struct ElasticShellDegenerate<:Material
    E::Float64
    nu::Float64
    t::Float64

    function ElasticShellDegenerate(prms::Dict{Symbol,Float64})
        return  ElasticShellDegenerate(;prms...)
    end

    function ElasticShellDegenerate(;E=NaN, nu=NaN, thickness=0.1)
        E>0.0 || error("Invalid value for E: $E")
        (0<=nu<0.5) || error("Invalid value for nu: $nu")
        thickness >0.0 || error("Invalid value for thickness: $thickness")
        this = new(E, nu, thickness)
        return this
    end
end

matching_elem_type(::ElasticShellDegenerate) = ShellDegenerate

# Type of corresponding state structure
ip_state_type(mat::ElasticShellDegenerate) = ElasticShellDegenerateIpState

function ip_state_vals(mat::ElasticShellDegenerate, ipd::ElasticShellDegenerateIpState)
    return OrderedDict{Symbol, Float64}()
end