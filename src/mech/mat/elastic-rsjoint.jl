# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticRSJoint

mutable struct ElasticRSJointState<:IpState
    env::ModelEnv
    σ ::Array{Float64,1}
    u ::Array{Float64,1}
    function ElasticRSJointState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ = zeros(env.ndim)
        this.u = zeros(env.ndim)
        return this
    end
end

mutable struct ElasticRSJoint<:MatParams
    ks::Float64
    kn::Float64

    function ElasticRSJoint(prms::Dict{Symbol,Float64})
        return  ElasticRSJoint(;prms...)
    end

    function ElasticRSJoint(;ks=NaN, kn=NaN)
        @check ks>=0
        @check kn>=0

        this = new(ks, kn)
        return this
    end
end


# Returns the element type that works with this material
matching_elem_type(::ElasticRSJoint) = MechRSJointElem

# Type of corresponding state structure
ip_state_type(::MechRSJointElem, ::ElasticRSJoint) = ElasticRSJointState

function calcD(matparams::ElasticRSJoint, state::ElasticRSJointState)
    ks = matparams.ks
    kn = matparams.kn
    if state.env.ndim==2
        return [  ks  0.0
                 0.0   kn ]
    else
        return  [  ks  0.0  0.0
                  0.0   kn  0.0
                  0.0  0.0   kn ]
    end
end


function update_state(matparams::ElasticRSJoint, state::ElasticRSJointState, Δu)
    D = calcD(matparams, state)
    Δσ = D*Δu

    state.u .+= Δu
    state.σ .+= Δσ
    return Δσ, success()
end

function ip_state_vals(matparams::ElasticRSJoint, state::ElasticRSJointState)
    return OrderedDict(
      :ur   => state.u[1] ,
      :tau  => state.σ[1] )
end
