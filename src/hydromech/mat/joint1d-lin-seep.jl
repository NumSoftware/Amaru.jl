# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export Joint1DLinSeep

mutable struct Joint1DLinSeepState<:IpState
    env::ModelEnv
    ndim::Int
    V::Float64     # fluid velocity
    D::Float64     # distance traveled by the fluid
    function Joint1DLinSeepState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.V = 0.0
        this.D = 0.0
        return this
    end
end

mutable struct Joint1DLinSeep<:Material
    k ::Float64    # specific permeability per meter

    function Joint1DLinSeep(; params...)
        names = (k="Permeability",)
        required = (:k,)
        @checkmissing params required names

        params  = (; params...)
        k       = params.k

        @check k>=0.0
        return new(k)
    end
end


# Element types that work with this material


# Type of corresponding state structure
ip_state_type(::Type{Joint1DLinSeep}) = Joint1DLinSeepState


function update_state!(mat::Joint1DLinSeep, state::Joint1DLinSeepState, ΔFw::Float64, Δt::Float64)
    k = mat.k
    state.V  = -k*ΔFw
    state.D  += state.V*Δt
    return state.V
end


function ip_state_vals(mat::Joint1DLinSeep, state::Joint1DLinSeepState)
    return OrderedDict(
      :vj => state.V)
end
