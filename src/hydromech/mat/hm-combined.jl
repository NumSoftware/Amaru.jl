
export HMCombined

mutable struct HMCombinedState{S1,S2}<:IpState
    env::ModelEnv
    hstate::S1 # hydro
    mstate::S2 # mech
    function HMCombinedState{S1,S2}(env=ModelEnv()) where {S1,S2}
        return new{S1,S2}(env, S1(env), S2(env))
    end
end


function Base.getproperty(state::HMCombinedState, s::Symbol)
    hstate = getfield(state, :hstate)
    mstate = getfield(state, :mstate)
    s==:hstate && return hstate
    s==:mstate && return mstate
    s==:env && return getfield(state, :env)
    fields1 = propertynames(hstate)
    fields2 = propertynames(mstate)
    s in fields1 && return getfield(hstate, s)
    s in fields2 && return getfield(mstate, s)
    error("type HMCombinedState has no field $s")
end


function Base.setproperty!(state::HMCombinedState, s::Symbol, value)
    hstate = getfield(state, :hstate)
    mstate = getfield(state, :mstate)
    fields1 = propertynames(hstate)
    fields2 = propertynames(mstate)
    
    s in fields1 && return setfield!(hstate, s, value)
    s in fields2 && return setfield!(hstate, s, value)
    error("type HMCombinedState has no field $s")
end


mutable struct HMCombined{M1,M2}<:Material
    tmat::M1 # hydro
    mmat::M2 # mech

    function HMCombined{M1,M2}(; params...) where {M1,M2}
        tmat = M1(;params...)
        mmat = M2(;params...)
        
        return new{M1,M2}(tmat, mmat)
    end
end

# Type of corresponding state structure
compat_state_type(::Type{HMCombined{M1,M2}}, ::Type{HMSolid}, env::ModelEnv) where {M1,M2} = HMCombinedState{compat_state_type(M1,SeepSolid,env), compat_state_type(M2,MechSolid,env)} 
compat_state_type(::Type{HMCombined{M1,M2}}, ::Type{HMJoint}, env::ModelEnv) where {M1,M2} = HMCombinedState{compat_state_type(M1,HydroJoint,env), compat_state_type(M2,MechJoint,env)} 

# compat_elem_types(::Type{HMCombined{M1,M2}}) where {M1,M2} = (HMSolid, HMJoint)


function Base.getproperty(mat::HMCombined, s::Symbol)
    tmat = getfield(mat, :tmat)
    mmat = getfield(mat, :mmat)
    s==:tmat && return tmat
    s==:mmat && return mmat
    fields1 = propertynames(tmat)
    fields2 = propertynames(mmat)
    s in fields1 && return getfield(tmat, s)
    s in fields2 && return getfield(mmat, s)
    error("type HMCombined has no field $s")
end


function calcD(mat::HMCombined, state::HMCombinedState)
    return calcD(mat.mmat, state.mstate)
end

function calcK(mat::HMCombined, state::HMCombinedState) # Hydraulic conductivity matrix
    return calcK(mat.tmat, state.hstate)
end

function update_state!(mat::HMCombined, state::HMCombinedState, Δε::Array{Float64,1}, Δut::Float64, G::Array{Float64,1}, Δt::Float64)
    QQ = update_state!(mat.tmat, state.hstate, Δut, G, Δt)
    Δσ, status = update_state!(mat.mmat, state.mstate, Δε)
    return Δσ, QQ, status
end

function ip_state_vals(mat::HMCombined, state::HMCombinedState)
    vals1 = ip_state_vals(mat.tmat, state.hstate)
    vals2 = ip_state_vals(mat.mmat, state.mstate)
    return merge(vals1, vals2)
end
