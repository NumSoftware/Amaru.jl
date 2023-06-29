# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticTipJoint

mutable struct ElasticTipJointState<:IpState
    env::ModelEnv
    f ::Float64
    w ::Float64
    function ElasticTipJointState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.f = 0.0
        this.w = 0.0
        return this
    end
end


mutable struct ElasticTipJoint<:Material
    k::Float64

    function ElasticTipJoint(prms::Dict{Symbol,Float64})
        return  ElasticTipJoint(;prms...)
    end

    function ElasticTipJoint(;k=NaN)
        @check k>=0
        this = new(k)
        return this
    end
end


# Returns the element type that works with this material
matching_elem_type(::ElasticTipJoint, shape::CellShape, ndim::Int) = MechTipJoint

# Type of corresponding state structure
ip_state_type(mat::ElasticTipJoint) = ElasticTipJointState

function calcD(mat::ElasticTipJoint, state::ElasticTipJointState)
    return mat.k
end

function stress_update(mat::ElasticTipJoint, state::ElasticTipJointState, Δw)
    Δf = mat.k*Δw
    state.f += Δf
    state.w += Δw
    return Δf, success()
end

function ip_state_vals(mat::ElasticTipJoint, state::ElasticTipJointState)
    return OrderedDict(
      :ur   => state.w ,
      :tau  => state.f )
end
