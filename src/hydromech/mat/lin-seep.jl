# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ConstPermeability

mutable struct ConstPermeabilityState<:IpState
    ctx::Context
    V::Array{Float64,1} # fluid velocity
    D::Array{Float64,1} # distance traveled by the fluid
    uw::Float64         # pore pressure
    function ConstPermeabilityState(ctx::Context)
        this = new(ctx)
        this.V  = zeros(ctx.ndim)
        this.D  = zeros(ctx.ndim)
        this.uw = 0.0
        return this
    end
end


mutable struct ConstPermeability<:Material
    k ::Float64 # specific permeability
    S ::Float64 # storativity coefficient

    function ConstPermeability(; params...)
        names = (k="Permeability", S="Storativity coefficient")
        required = (:k, )
        @checkmissing params required names

        default = (S=0.0,)
        params  = merge(default, params)

        params  = (; params...)
        k       = params.k
        S       = params.S

        @check k>0.0
        @check S>=0.0

        return new(k, S)
    end
end


# Type of corresponding state structure
compat_state_type(::Type{ConstPermeability}, ::Type{SeepSolid}, ctx::Context) = ConstPermeabilityState

# Element types that work with this material
# compat_elem_types(::Type{ConstPermeability}) = (SeepSolid,)


function calcK(mat::ConstPermeability, state::ConstPermeabilityState) # Hydraulic conductivity matrix
    if state.ctx.ndim==2
        return mat.k*eye(2)
    else
        return mat.k*eye(3)
    end
end


function update_state!(mat::ConstPermeability, state::ConstPermeabilityState, Δuw::Float64, G::Array{Float64,1}, Δt::Float64)
    K = calcK(mat, state)
    state.V   = -K*G
    state.D  += state.V*Δt
    state.uw += Δuw
    return state.V
end


function ip_state_vals(mat::ConstPermeability, state::ConstPermeabilityState)
    D = OrderedDict{Symbol, Float64}()
    D[:vx] = state.V[1]
    D[:vy] = state.V[2]
    if state.ctx.ndim==3
        D[:vz] = state.V[3]
    end

    return D
end
