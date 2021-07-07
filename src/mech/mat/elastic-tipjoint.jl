# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticTipJoint

mutable struct ElasticTipJointIpState<:IpState
    env::ModelEnv
    f ::Float64
    w ::Float64
    function ElasticTipJointIpState(env::ModelEnv=ModelEnv())
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
matching_elem_type(::ElasticTipJoint) = MechTipJoint

# Type of corresponding state structure
ip_state_type(mat::ElasticTipJoint) = ElasticTipJointIpState

function calcD(mat::ElasticTipJoint, ipd::ElasticTipJointIpState)
    return mat.k
end

function stress_update(mat::ElasticTipJoint, ipd::ElasticTipJointIpState, Δw)
    Δf = mat.k*Δw
    ipd.f += Δf
    ipd.w += Δw
    return Δf, success()
end

function ip_state_vals(mat::ElasticTipJoint, ipd::ElasticTipJointIpState)
    return OrderedDict(
      :ur   => ipd.w ,
      :tau  => ipd.f )
end
