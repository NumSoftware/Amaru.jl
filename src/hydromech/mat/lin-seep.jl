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


mutable struct LinSeep<:MatParams
    k ::Float64 # specific permeability
    S ::Float64 # storativity coefficient

    function LinSeep(;k=NaN, S=0.0)
        @check k>0.0
        @check S>=0.0

        this = new(k, S)
        return this
    end
end



# Type of corresponding state structure
ip_state_type(::SeepSolidElem, ::LinSeep) = LinSeepState


function calcK(matparams::LinSeep, state::LinSeepState) # Hydraulic conductivity matrix
    if state.env.ndim==2
        return matparams.k*eye(2)
    else
        return matparams.k*eye(3)
    end
end


function update_state!(matparams::LinSeep, state::LinSeepState, Δuw::Float64, G::Array{Float64,1}, Δt::Float64)
    K = calcK(matparams, state)
    state.V   = -K*G
    state.D  += state.V*Δt
    state.uw += Δuw
    return state.V
end


function ip_state_vals(matparams::LinSeep, state::LinSeepState)
    D = OrderedDict{Symbol, Float64}()
    D[:vx] = state.V[1]
    D[:vy] = state.V[2]
    if state.env.ndim==3
        D[:vz] = state.V[3]
    end

    return D
end
