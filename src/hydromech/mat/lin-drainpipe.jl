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

mutable struct LinDrainPipe<:MatParams
    k ::Float64     # specific permeability

    function LinDrainPipe(prms::Dict{Symbol,Float64})
        return  LinDrainPipe(;prms...)
    end

    function LinDrainPipe(;k=NaN)
        k<=0.0 && error("Invalid value for k: $k")
        return new(k)
    end
end

matching_elem_type(::LinDrainPipe) = DrainPipeElem
#matching_elem_type_if_embedded(::LinDrainPipe) = SeepEmbRod

# Type of corresponding state structure
ip_state_type(::DrainPipeElem, ::LinDrainPipe) = LinDrainPipeState


function update_state!(matparams::LinDrainPipe, state::LinDrainPipeState, Δuw::Float64, G::Float64, Δt::Float64)
    k = matparams.k
    state.V  = -k*G
    state.D  += state.V*Δt
    state.uw += Δuw
    return state.V
end


function ip_state_vals(matparams::LinDrainPipe, state::LinDrainPipeState)
    return OrderedDict(
      :va => state.V,
      :uwa => state.uw)
    #   :Fa => state.V*matparams.A,
    #   :A  => matparams.A )
end


