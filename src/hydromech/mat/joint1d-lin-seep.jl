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

mutable struct Joint1DLinSeep<:MatParams
    k ::Float64    # specific permeability per meter

    function Joint1DLinSeep(prms::Dict{Symbol,Float64})
        return  Joint1DLinSeep(;prms...)
    end

    function Joint1DLinSeep(;k=NaN)
        @check k>=0.0
        return new(k)
    end
end


# Returns the element type that works with this material


# Type of corresponding state structure
ip_state_type(::SeepJoint1DElem, ::Joint1DLinSeep) = Joint1DLinSeepState


function update_state!(matparams::Joint1DLinSeep, state::Joint1DLinSeepState, ΔFw::Float64, Δt::Float64)
    k = matparams.k
    state.V  = -k*ΔFw
    state.D  += state.V*Δt
    return state.V
end


function ip_state_vals(matparams::Joint1DLinSeep, state::Joint1DLinSeepState)
    return OrderedDict(
      :vj => state.V)
end
