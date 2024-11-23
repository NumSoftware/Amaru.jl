# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticPlateMZC

mutable struct ElasticPlateMZCState<:IpState
    ctx::Context
    function ElasticPlateMZCState(ctx::Context)
        return new(ctx)
    end
end

mutable struct ElasticPlateMZC<:Material
    E::Float64
    ν::Float64
    thick::Float64
    ρ::Float64

    function ElasticPlateMZC(prms::Dict{Symbol,Float64})
        return  ElasticPlateMZC(;prms...)
    end

    function ElasticPlateMZC(;E=NaN, nu=NaN, thick=NaN, ρ=0.0)
        E>0.0 || error("Invalid value for E: $E")
        (0<=nu<0.5) || error("Invalid value for nu: $nu")
        thick>0.0 || error("Invalid value for thick: $thick")
        this = new(E, nu, thick, ρ)
        return this
    end
end



# Type of corresponding state structure
compat_state_type(::Type{ElasticPlateMZC}) = ElasticPlateMZCState


function ip_state_vals(mat::ElasticPlateMZC, state::ElasticPlateMZCState)
    return OrderedDict{Symbol, Float64}()
end
