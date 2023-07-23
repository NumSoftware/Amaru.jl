# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LinDrainPipe

mutable struct LinDrainPipeState<:IpState
    env::ModelEnv
    V::Float64       # fluid velocity
    D::Float64       # distance traveled by the fluid
    uw::Float64      # pore pressure
    function LinDrainPipeState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.V  = 0.0
        this.D  = 0.0
        this.uw = 0.0
        return this
    end
end

mutable struct LinDrainPipe<:Material
    k ::Float64     # specific permeability

    function LinDrainPipe(; params...)
        names = (k="Permeability",)
        required = (:k,)
        @checkmissing params required names

        params  = (; params...)
        k       = params.k
        @check k>0.0

        return new(k)
    end
end


# Type of corresponding state structure
compat_state_type(::Type{LinDrainPipe}) = LinDrainPipeState

# Element types that work with this material
compat_elem_types(::Type{LinDrainPipe}) = (DrainPipe,)


function update_state!(mat::LinDrainPipe, state::LinDrainPipeState, Δuw::Float64, G::Float64, Δt::Float64)
    k = mat.k
    state.V  = -k*G
    state.D  += state.V*Δt
    state.uw += Δuw
    return state.V
end


function ip_state_vals(mat::LinDrainPipe, state::LinDrainPipeState)
    return OrderedDict(
      :va => state.V,
      :uwa => state.uw)
    #   :Fa => state.V*mat.A,
    #   :A  => mat.A )
end


