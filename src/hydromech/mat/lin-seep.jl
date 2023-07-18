# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LinSeep

mutable struct LinSeepState<:IpState
    env::ModelEnv
    V::Array{Float64,1} # fluid velocity
    D::Array{Float64,1} # distance traveled by the fluid
    uw::Float64         # pore pressure
    function LinSeepState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.V  = zeros(env.ndim)
        this.D  = zeros(env.ndim)
        this.uw = 0.0
        return this
    end
end


mutable struct LinSeep<:Material
    k ::Float64 # specific permeability
    S ::Float64 # storativity coefficient

    function LinSeep(; params...)
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
ip_state_type(::SeepSolid, ::LinSeep) = LinSeepState


function calcK(mat::LinSeep, state::LinSeepState) # Hydraulic conductivity matrix
    if state.env.ndim==2
        return mat.k*eye(2)
    else
        return mat.k*eye(3)
    end
end


function update_state!(mat::LinSeep, state::LinSeepState, Δuw::Float64, G::Array{Float64,1}, Δt::Float64)
    K = calcK(mat, state)
    state.V   = -K*G
    state.D  += state.V*Δt
    state.uw += Δuw
    return state.V
end


function ip_state_vals(mat::LinSeep, state::LinSeepState)
    D = OrderedDict{Symbol, Float64}()
    D[:vx] = state.V[1]
    D[:vy] = state.V[2]
    if state.env.ndim==3
        D[:vz] = state.V[3]
    end

    return D
end
