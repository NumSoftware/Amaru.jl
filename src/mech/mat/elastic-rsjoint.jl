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
    p ::Float64    # section perimeter

    function ElasticRSJoint(prms::Dict{Symbol,Float64})
        return  ElasticRSJoint(;prms...)
    end

    function ElasticRSJoint(;ks=NaN, kn=NaN, p=NaN, A=NaN, dm=NaN)
        # A : section area
        # dm: section diameter
        # p : section perimeter
        @check ks>=0
        @check kn>=0
        @check (p>0 || A>0 || dm>0)

        if isnan(p)
            if A>0
                p = 2.0*(A*pi)^0.5
            else
                p = pi*dm
            end
        end
        @assert p>0

        this = new(ks, kn, p)
        return this
    end
end

ElasticJoint1D = ElasticRSJoint #! deprecated
export ElasticJoint1D


# Returns the element type that works with this material
matching_elem_type(::ElasticRSJoint, shape::CellShape, ndim::Int) = MechRodSolidJoint

# Type of corresponding state structure
ip_state_type(mat::ElasticRSJoint) = ElasticRSJointState

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


function stress_update(mat::ElasticRSJoint, state::ElasticRSJointState, Δu)
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
