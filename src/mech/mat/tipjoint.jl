# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export TipJoint

mutable struct TipJointState<:IpState
    env::ModelEnv
    f ::Float64
    w ::Float64
    function TipJointState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.f = 0.0
        this.w = 0.0
        return this
    end
end


mutable struct TipJoint<:Material
    k::Float64

    function TipJoint(prms::Dict{Symbol,Float64})
        return  TipJoint(;prms...)
    end

    function TipJoint(;k=NaN)
        @check k>=0
        this = new(k)
        return this
    end
end


# Element types that work with this material
compat_elem_types(::Type{TipJoint}) = (MechTipJoint,)

# Type of corresponding state structure
compat_state_type(::Type{TipJoint}) = TipJointState


function calcD(mat::TipJoint, state::TipJointState)
    if state.w>0.0
        return mat.k
    else
        return 0.0
    end
end


function update_state!(mat::TipJoint, state::TipJointState, Δw)
    fini = state.f
    ftr  = fini + mat.k*Δw
    
    if ftr>0.0
        f = ftr
    else
        f = 0.0
    end
    Δf = f - fini
    state.f += Δf
    state.w += Δw
    return Δf, success()
end


function ip_state_vals(mat::TipJoint, state::TipJointState)
    return OrderedDict(
      :ur   => state.w ,
      :tau  => state.f )
end
