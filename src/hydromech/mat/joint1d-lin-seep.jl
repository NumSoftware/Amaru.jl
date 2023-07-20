# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export Joint1DConstPermeability

mutable struct Joint1DConstPermeabilityState<:IpState
    env::ModelEnv
    ndim::Int
    V::Float64     # fluid velocity
    D::Float64     # distance traveled by the fluid
    function Joint1DConstPermeabilityState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.V = 0.0
        this.D = 0.0
        return this
    end
end

mutable struct Joint1DConstPermeability<:Material
    k ::Float64    # specific permeability per meter

    function Joint1DConstPermeability(; params...)
        names = (k="Permeability",)
        required = (:k,)
        @checkmissing params required names

        params  = (; params...)
        k       = params.k

        @check k>=0.0
        return new(k)
    end
end


# Type of corresponding state structure
ip_state_type(::Type{Joint1DConstPermeability}) = Joint1DConstPermeabilityState

# Element types that work with this material
matching_elem_types(::Type{Joint1DConstPermeability}) = (SeepJoint1D,)


function update_state!(mat::Joint1DConstPermeability, state::Joint1DConstPermeabilityState, ΔFw::Float64, Δt::Float64)
    k = mat.k
    state.V  = -k*ΔFw
    state.D  += state.V*Δt
    return state.V
end


function ip_state_vals(mat::Joint1DConstPermeability, state::Joint1DConstPermeabilityState)
    return OrderedDict(
      :vj => state.V)
end
