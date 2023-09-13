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

mutable struct ElasticRSJoint<:Material
    ks::Float64
    kn::Float64

    function ElasticRSJoint(prms::Dict{Symbol,Float64})
        return  ElasticRSJoint(;prms...)
    end

    function ElasticRSJoint(; params...)
        names = (kn="Normal stiffness", ks="Shear stiffness")
        required = keys(names)
        @checkmissing params required names
        
        params = (; params...)
        kn     = params.kn
        ks     = params.ks

        @check ks>=0
        @check kn>=0

        return new(ks, kn)
    end
end


# Type of corresponding state structure
compat_state_type(::Type{ElasticRSJoint}, ::Type{MechRSJoint}, env::ModelEnv) = ElasticRSJointState

# Element types that work with this material
# compat_elem_types(::Type{ElasticRSJoint}) = (MechRSJoint,)


function calcD(mat::ElasticRSJoint, state::ElasticRSJointState)
    ks = mat.ks
    kn = mat.kn
    if state.env.ndim==2
        return [  ks  0.0
                 0.0   kn ]
    else
        return  [  ks  0.0  0.0
                  0.0   kn  0.0
                  0.0  0.0   kn ]
    end
end


function update_state!(mat::ElasticRSJoint, state::ElasticRSJointState, Δu)
    D = calcD(mat, state)
    Δσ = D*Δu

    state.u .+= Δu
    state.σ .+= Δσ
    return Δσ, success()
end


function ip_state_vals(mat::ElasticRSJoint, state::ElasticRSJointState)
    return OrderedDict(
      :ur   => state.u[1] ,
      :tau  => state.σ[1] )
end
