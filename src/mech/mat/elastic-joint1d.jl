# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticLineJoint

mutable struct ElasticLineJointIpState<:IpState
    env::ModelEnv
    σ ::Array{Float64,1}
    u ::Array{Float64,1}
    function ElasticLineJointIpState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ = zeros(env.ndim)
        this.u = zeros(env.ndim)
        return this
    end
end


mutable struct ElasticLineJoint<:Material
    ks::Float64
    kn::Float64
    p ::Float64    # section perimeter

    function ElasticLineJoint(prms::Dict{Symbol,Float64})
        return  ElasticLineJoint(;prms...)
    end

    function ElasticLineJoint(;ks=NaN, kn=NaN, p=NaN, A=NaN, dm=NaN)
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

ElasticJoint1D = ElasticLineJoint #! deprecated
export ElasticJoint1D


# Returns the element type that works with this material
matching_elem_type(::ElasticLineJoint) = MechLineJoint

# Type of corresponding state structure
ip_state_type(mat::ElasticLineJoint) = ElasticLineJointIpState

function calcD(mat::ElasticLineJoint, ipd::ElasticLineJointIpState)
    ks = mat.ks
    kn = mat.kn
    if ipd.env.ndim==2
        return [  ks  0.0
                 0.0   kn ]
    else
        return  [  ks  0.0  0.0
                  0.0   kn  0.0
                  0.0  0.0   kn ]
    end
end


function stress_update(mat::ElasticLineJoint, ipd::ElasticLineJointIpState, Δu)
    D = calcD(mat, ipd)
    Δσ = D*Δu

    ipd.u .+= Δu
    ipd.σ .+= Δσ
    return Δσ, success()
end

function ip_state_vals(mat::ElasticLineJoint, ipd::ElasticLineJointIpState)
    return OrderedDict(
      :ur   => ipd.u[1] ,
      :tau  => ipd.σ[1] )
end
