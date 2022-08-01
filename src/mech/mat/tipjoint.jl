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


# Returns the element type that works with this material
matching_elem_type(::TipJoint) = MechTipJoint

# Type of corresponding state structure
ip_state_type(mat::TipJoint) = TipJointState

function calcD(mat::TipJoint, ipd::TipJointState)
    if ipd.w>0.0
        return mat.k
    else
        return 0.0
    end
end

function stress_update(mat::TipJoint, ipd::TipJointState, Δw)
    fini = ipd.f
    ftr  = fini + mat.k*Δw
    
    if ftr>0.0
        f = ftr
    else
        f = 0.0
    end
    Δf = f - fini
    ipd.f += Δf
    ipd.w += Δw
    return Δf, success()
end

function ip_state_vals(mat::TipJoint, ipd::TipJointState)
    return OrderedDict(
      :ur   => ipd.w ,
      :tau  => ipd.f )
end
