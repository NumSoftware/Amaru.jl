# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LinPress

mutable struct LinPressState<:IpState
    env::ModelEnv
    V  ::Array{Float64,1} # fluid velocity
    up ::Float64          # pressure
    function LinPressState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.V  = zeros(env.ndim)
        this.up = 0.0
        return this
    end
end


mutable struct LinPress<:Material
    μ::Float64 # viscocity
    c::Float64 # sound speed
    ρ::Float64 # densidade

    function LinPress(prms::Dict{Symbol,Float64})
        return  LinPress(;prms...)
    end

    function LinPress(;mu=NaN, rho=NaN, c=0.0)
        @check mu>0
        @check rho>0
        @check c>0

        this    = new(k, gammaw, S)
        return this
    end
end

# Returns the element type that works with this material model
matching_elem_type(::LinPress) = PressSolid

# Type of corresponding state structure
ip_state_type(mat::LinPress) = LinPressState



function update_state!(mat::LinPress, state::LinPressState, Δuw::Float64, G::Array{Float64,1}, Δt::Float64)
    K = calcK(mat, state)
    state.V   = -K*G
    state.D  += state.V*Δt
    state.up += Δuw
    return state.V
end


function ip_state_vals(mat::LinPress, state::LinPressState)
    D = OrderedDict{Symbol, Float64}()
    D[:up] = state.up
    D[:vx] = state.V[1]
    D[:vy] = state.V[2]
    if state.env.ndim==3
        D[:vz] = state.V[3]
    end

    return D
end
