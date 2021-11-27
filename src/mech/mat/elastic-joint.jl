# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticJoint

mutable struct JointIpState<:IpState
    env::ModelEnv
    σ   ::Array{Float64,1}
    w   ::Array{Float64,1}
    h   ::Float64
    function JointIpState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ   = zeros(3)
        this.w = zeros(3)
        this.h = 0.0
        return this
    end
end

mutable struct ElasticJoint<:Material
    E::Float64 # Young modulus from bulk material
    ν::Float64 # Poisson ration from bulk material
    ζ::Float64 # elastic displacement scale factor (formerly α)

    function ElasticJoint(prms::Dict{Symbol,Float64})
        return  ElasticJoint(;prms...)
    end

    function ElasticJoint(;E=NaN, nu=NaN, zeta=1.0)
        E>0.0       || error("Invalid value for E: $E")
        0<=nu<0.5   || error("Invalid value for nu: $nu")
        zeta>0      || error("Invalid value for zeta: $zeta")

        this = new(E, nu, zeta)
        return this
    end
end

# Returns the element type that works with this material model
@static if @isdefined MechJoint
    matching_elem_type(::ElasticJoint) = MechJoint
end

# Type of corresponding state structure
ip_state_type(mat::ElasticJoint) = JointIpState


function mountD(mat::ElasticJoint, ipd::JointIpState)
    ndim = ipd.env.ndim
    G  = mat.E/(1.0+mat.ν)/2.0
    kn = mat.E*mat.ζ/ipd.h
    ks =     G*mat.ζ/ipd.h
    if ndim==2
        return [  kn  0.0
                 0.0   ks ]
    else
        return  [  kn  0.0  0.0
                  0.0   ks  0.0
                  0.0  0.0   ks ]
    end
end

function stress_update(mat::ElasticJoint, ipd::JointIpState, Δu)
    ndim = ipd.env.ndim
    D  = mountD(mat, ipd)
    Δσ = D*Δu

    ipd.w[1:ndim] += Δu
    ipd.σ[1:ndim] += Δσ
    return Δσ, success()
end


function ip_state_vals(mat::ElasticJoint, ipd::JointIpState)
    ndim = ipd.env.ndim
    if ndim == 3
       return Dict(
          :wj1  => ipd.w[1],
          :wj2  => ipd.w[2],
          :wj3  => ipd.w[3],
          :sj1  => ipd.σ[1],
          :sj2  => ipd.σ[2],
          :sj3  => ipd.σ[3],
          )
    else
        return Dict(
          :wj1  => ipd.w[1],
          :wj2  => ipd.w[2],
          :sj1  => ipd.σ[1],
          :sj2  => ipd.σ[2],
          )
    end
end

function output_keys(mat::ElasticJoint)
    return Symbol[:wj1, :sj1]
end