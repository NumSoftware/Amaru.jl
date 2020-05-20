# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticShellQuad4node

mutable struct ElasticShellQuad4nodeIpState<:IpState
    env::ModelEnv
    function ElasticShellQuad4nodeIpState(env::ModelEnv=ModelEnv())
        return new(env)
    end
end

mutable struct ElasticShellQuad4node<:Material
    E::Float64
    nu::Float64
    thick::Float64

    function ElasticShellQuad4node(prms::Dict{Symbol,Float64})
        return  ElasticShellQuad4node(;prms...)
    end

    function ElasticShellQuad4node(;E=NaN, nu=NaN, thick=0.1)
        E>0.0 || error("Invalid value for E: $E")
        (0<=nu<0.5) || error("Invalid value for nu: $nu")
        thick >0.0 || error("Invalid value for thick: $thick")
        this = new(E, nu, thick)
        return this
    end
end

matching_elem_type(::ElasticShellQuad4node) = ShellQuad4node

# Type of corresponding state structure
ip_state_type(mat::ElasticShellQuad4node) = ElasticShellQuad4nodeIpState

function ip_state_vals(mat::ElasticShellQuad4node, ipd::ElasticShellQuad4nodeIpState)
    return OrderedDict{Symbol, Float64}()
end
