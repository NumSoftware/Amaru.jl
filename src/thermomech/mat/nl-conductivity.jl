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


mutable struct NLConductivity<:Material
    #   T   C
    #   0   55
    #   800  28
    #   1200 28
    k::Float64
    cv::Float64
    k_table::Array{Float64,2}
    cv_table::Array{Float64,2}

    function NLConductivity(; args...)
        args = checkargs(args, arg_rules(NLConductivity))

        if args.k isa Array
            k = 0.0
            k_table = args.k
        else
            k = args.k
            k_table = zeros(0,0)
        end

        if args.cv isa Array
            cv = 0.0
            cv_table = args.cv
        else
            cv = args.cv
            cv_table = zeros(0,0)
        end

        return new(k, cv, k_table, cv_table)
    end
end


arg_rules(::Type{NLConductivity}) =
[
    @arginfo k 1==1 "Conductivity"
    @arginfo cv=0 1==1 "Specific heat"
]

# Type of corresponding state structure
compat_state_type(::Type{NLConductivity}, ::Type{ThermoSolid}, env::ModelEnv) = NLConductivityState
compat_state_type(::Type{NLConductivity}, ::Type{ThermoShell}, env::ModelEnv) = NLConductivityState

# Element types that work with this material
# compat_elem_types(::Type{NLConductivity}) = (ThermoSolid,)


function calc_cv(mat::NLConductivity, ut::Float64) # Specific heat
    length(mat.cv_table)==0 && return mat.cv

    T  = mat.cv_table[:,1]
    Cv = mat.cv_table[:,2]

    i = searchsortedfirst(T, ut)
    if i==1
        cv = Cv[1]
    elseif i>length(T)
        cv = Cv[end]
    else
        cv = Cv[i-1] + (ut-T[i-1]) * (Cv[i]-Cv[i-1])/(T[i]-T[i-1])
    end

    return cv
end

function calcK(mat::NLConductivity, ut::Float64) # Thermal conductivity matrix
    length(mat.k_table)==0 && return mat.k

    T  = mat.k_table[:,1]
    C = mat.k_table[:,2]

    i = searchsortedfirst(T, ut)
    if i==1
        k = C[1]
    elseif i>length(T)
        k = C[end]
    else
        k = C[i-1] + (ut-T[i-1]) * (C[i]-C[i-1])/(T[i]-T[i-1])
    end

    return k
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
