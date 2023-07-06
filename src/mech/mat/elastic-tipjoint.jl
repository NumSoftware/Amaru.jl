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


mutable struct ElasticTipJoint<:MatParams
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
matching_elem_type(::ElasticTipJoint) = MechTipJointElem

# Type of corresponding state structure
ip_state_type(matparams::ElasticTipJoint) = ElasticTipJointState

function calcD(matparams::ElasticTipJoint, state::ElasticTipJointState)
    return matparams.k
end

function update_state(matparams::ElasticTipJoint, state::ElasticTipJointState, Δw)
    Δf = matparams.k*Δw
    state.f += Δf
    state.w += Δw
    return Δf, success()
end

function ip_state_vals(matparams::ElasticTipJoint, state::ElasticTipJointState)
    return OrderedDict(
      :ur   => state.w ,
      :tau  => state.f )
end
