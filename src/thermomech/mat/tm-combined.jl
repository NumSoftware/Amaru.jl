
export TMCombined

mutable struct TMCombinedState{S1,S2}<:IpState
    env::ModelEnv
    tstate::S1 # thermo
    mstate::S2 # mech
    function TMCombinedState{S1,S2}(env=ModelEnv()) where {S1,S2}
        return new{S1,S2}(env, S1(env), S2(env))
    end
end


function Base.getproperty(state::TMCombinedState, s::Symbol)
    tstate = getfield(state, :tstate)
    mstate = getfield(state, :mstate)
    s==:tstate && return tstate
    s==:mstate && return mstate
    s==:env && return getfield(state, :env)
    fields1 = propertynames(tstate)
    fields2 = propertynames(mstate)
    s in fields1 && return getfield(tstate, s)
    s in fields2 && return getfield(mstate, s)
    error("type TMCombinedState has no field $s")
end


TMCombined_params = [
    FunInfo(:TMCombined, "A material model for combined thermal and mechanical effects"),
    KwArgInfo(:alpha, "Thermal expansion coefficient", cond=:(alpha>=0))
]
@doc docstring(TMCombined_params) TMCombined

mutable struct TMCombined{M1,M2}<:Material
    tmat::M1 # thermo
    mmat::M2 # mech
    α::Float64 # thermal expansion coefficient
    α_table::Array{Float64,2}

    function TMCombined{M1,M2}(; args...) where {M1,M2}
        tmat = M1(;args...)
        mmat = M2(;args...)

        args = checkargs(args, TMCombined_params)

        if args.alpha isa Array
            alpha = 0.0
            alpha_table = args.alpha
        else
            alpha = args.alpha
            alpha_table = zeros(0,0)
        end

        return new{M1,M2}(tmat, mmat, alpha, alpha_table)
    end
end


function calc_α(mat::TMCombined{M1,M2}, ut::Float64) where {M1,M2} # thermal expansion coefficient  1/K or 1/°C
    length(mat.α_table)==0 && return mat.α

    T = elem.props.α_table[:,1]
    A = elem.props.α_table[:,2]

    i = searchsortedfirst(T, ut)
    if i==1
        α = A[1]
    elseif i>length(T)
        α = A[end]
    else
        α = A[i-1] + (ut-T[i-1]) * (A[i]-A[i-1])/(T[i]-T[i-1])
    end

    return α
end


# Type of corresponding state structure
compat_state_type(::Type{TMCombined{M1,M2}}, ::Type{TMSolid}, env::ModelEnv) where {M1,M2} = TMCombinedState{compat_state_type(M1,ThermoSolid,env), compat_state_type(M2,MechSolid,env)} 
compat_state_type(::Type{TMCombined{M1,M2}}, ::Type{TMShell}, env::ModelEnv) where {M1,M2} = TMCombinedState{compat_state_type(M1,ThermoShell,env), compat_state_type(M2,MechShell,env)} 


function Base.getproperty(mat::TMCombined, s::Symbol)
    tmat = getfield(mat, :tmat)
    mmat = getfield(mat, :mmat)
    s==:tmat && return tmat
    s==:mmat && return mmat
    s in propertynames(tmat) && return getfield(tmat, s)
    s in propertynames(mmat) && return getfield(mmat, s)
    s in propertynames(mat) && return getfield(mat, s)
    error("type TMCombined has no field $s")
end


function calcD(mat::TMCombined, state::TMCombinedState)
    return calcD(mat.mmat, state.mstate)
end

function calc_cv(mat::TMCombined, ut::Float64) # Hydraulic conductivity matrix
    return calc_cv(mat.tmat, ut)
end

function calcK(mat::TMCombined, state::TMCombinedState) # Hydraulic conductivity matrix
    return calcK(mat.tmat, state.tstate)
end

function update_state!(mat::TMCombined, state::TMCombinedState, Δε::Array{Float64,1}, Δut::Float64, G::Array{Float64,1}, Δt::Float64)
    QQ = update_state!(mat.tmat, state.tstate, Δut, G, Δt)
    Δσ, status = update_state!(mat.mmat, state.mstate, Δε)
    return Δσ, QQ, status
end

function ip_state_vals(mat::TMCombined, state::TMCombinedState)
    vals1 = ip_state_vals(mat.tmat, state.tstate)
    vals2 = ip_state_vals(mat.mmat, state.mstate)
    return merge(vals1, vals2)
end
