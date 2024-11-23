# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticShellDegenerated

mutable struct ElasticShellDegeneratedState<:IpState
    ctx::Context
    σ::Array{Float64,1}
    function ElasticShellDegeneratedState(ctx::Context)
        σ = zeros(5)
        return new(ctx, σ)
    end
end

mutable struct ElasticShellDegenerated<:Material
    E::Float64
    ν::Float64
    t::Float64

    function ElasticShellDegenerated(prms::Dict{Symbol,Float64})
        return  ElasticShellDegenerated(;prms...)
    end

    function ElasticShellDegenerated(;E=NaN, nu=NaN, thickness=0.1)
        E>0.0 || error("Invalid value for E: $E")
        (0<=nu<0.5) || error("Invalid value for nu: $nu")
        thickness >0.0 || error("Invalid value for thickness: $thickness")
        this = new(E, nu, thickness)
        return this
    end
end



# Type of corresponding state structure
compat_state_type(::Type{ElasticShellDegenerated}) = ElasticShellDegeneratedState


function ip_state_vals(mat::ElasticShellDegenerated, state::ElasticShellDegeneratedState)
    return OrderedDict{Symbol, Float64}()
end