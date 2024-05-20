# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export NLConductivity

mutable struct NLConductivityState<:IpState
    env::ModelEnv
    ut::Float64
    Q::Array{Float64,1} # thermal flow
    D::Array{Float64,1}
    function NLConductivityState(env::ModelEnv)
        this = new(env)
        this.ut = 0.0
        this.Q  = zeros(env.ndim)
        this.D  = zeros(env.ndim)
        return this
    end
end


NLConductivity_params = [
    FunInfo(:NLConductivity, "Nonlinear thermal conductivity material model"),
    KwArgInfo(:k, "Conductivity"),
    KwArgInfo(:cv, "Specific heat", 0.0)
]
@doc docstring(NLConductivity_params) NLConductivity

mutable struct NLConductivity<:Material
    #   T   C
    #   0   55
    #   800  28
    #   1200 28
    k::Float64
    cv::Float64
    k_fun::PathFunction
    cv_fun::PathFunction

    function NLConductivity(; args...)
        args = checkargs(args, NLConductivity_params)
        this = new()

        if args.k isa PathFunction
            this.k = 0.0
            this.k_fun = args.k
        else
            this.k = args.k
        end

        if args.cv isa PathFunction
            this.cv = 0.0
            this.cv_fun = args.cv
        else
            this.cv = args.cv
        end

        return this
    end
end


# Type of corresponding state structure
compat_state_type(::Type{NLConductivity}, ::Type{ThermoSolid}, env::ModelEnv) = NLConductivityState
compat_state_type(::Type{NLConductivity}, ::Type{ThermoShell}, env::ModelEnv) = NLConductivityState


function calc_cv(mat::NLConductivity, ut::Float64) # Specific heat
    isdefined(mat, :cv_fun) || return mat.cv

    return mat.cv_fun(ut)
end

function calcK(mat::NLConductivity, state::NLConductivityState) # Thermal conductivity
    isdefined(mat, :k_fun) || return mat.k

    return mat.k_fun(state.ut)
end


function update_state!(mat::NLConductivity, state::NLConductivityState, Δut::Float64, G::Array{Float64,1}, Δt::Float64)
    K = calcK(mat, state)
    state.Q   = -K*G
    state.D  += state.Q*Δt
    state.ut += Δut
    return state.Q
end


function ip_state_vals(mat::NLConductivity, state::NLConductivityState)
    D = OrderedDict{Symbol, Float64}()
    #D[:qx] = state.Q[1]
    #D[:qy] = state.Q[2]
    #if state.env.ndim==3
        #D[:qz] = state.Q[3]
    #end
    return D
end
