
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


mutable struct TMCombined{M1,M2}<:Material
    tmat::M1 # thermo
    mmat::M2 # mech

    function TMCombined{M1,M2}(; params...) where {M1,M2}
        tmat = M1(;params...)
        mmat = M2(;params...)
        
        return new{M1,M2}(tmat, mmat)
    end
end

# Type of corresponding state structure
compat_state_type(::Type{TMCombined{M1,M2}}, ::Type{TMSolid}, env::ModelEnv) where {M1,M2} = TMCombinedState{compat_state_type(M1,ThermoSolid,env), compat_state_type(M2,MechSolid,env)} 
compat_state_type(::Type{TMCombined{M1,M2}}, ::Type{TMShell}, env::ModelEnv) where {M1,M2} = TMCombinedState{compat_state_type(M1,ThermoShell,env), compat_state_type(M2,MechShell,env)} 

# compat_elem_types(::Type{TMCombined{M1,M2}}) where {M1,M2} = (TMSolid, TMShell)


function Base.getproperty(mat::TMCombined, s::Symbol)
    tmat = getfield(mat, :tmat)
    mmat = getfield(mat, :mmat)
    s==:tmat && return tmat
    s==:mmat && return mmat
    fields1 = propertynames(tmat)
    fields2 = propertynames(mmat)
    s in fields1 && return getfield(tmat, s)
    s in fields2 && return getfield(mmat, s)
    error("type TMCombined has no field $s")
end


function calcD(mat::TMCombined, state::TMCombinedState, stressmodel::String=state.env.ana.stressmodel)
    return calcD(mat.mmat, state.mstate, stressmodel)
end

function calc_cv(mat::TMCombined, ut::Float64) # Hydraulic conductivity matrix
    return calc_cv(mat.tmat, ut)
end

function calcK(mat::TMCombined, state::TMCombinedState) # Hydraulic conductivity matrix
    return calcK(mat.tmat, state.tstate)
end

function update_state!(mat::TMCombined, state::TMCombinedState, Δε::Array{Float64,1}, Δut::Float64, G::Array{Float64,1}, Δt::Float64, stressmodel::String=state.env.ana.stressmodel)
    QQ = update_state!(mat.tmat, state.tstate, Δut, G, Δt)
    Δσ, status = update_state!(mat.mmat, state.mstate, Δε, stressmodel)
    return Δσ, QQ, status
end

function ip_state_vals(mat::TMCombined, state::TMCombinedState, stressmodel::String=state.env.ana.stressmodel)
    vals1 = ip_state_vals(mat.tmat, state.tstate)
    vals2 = ip_state_vals(mat.mmat, state.mstate, stressmodel)
    return merge(vals1, vals2)
end
