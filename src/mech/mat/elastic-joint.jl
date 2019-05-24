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
    ζ::Float64 # elastic displacement scale factor

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

# Create a new instance of Ip data
new_ip_state(mat::ElasticJoint, env::ModelEnv) = JointIpState(env)


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
    return Δσ
end

function ip_state_vals(mat::ElasticJoint, ipd::JointIpState)
    ndim = ipd.env.ndim
    if ndim == 2
        return OrderedDict(
          :w1  => ipd.w[1] ,
          :w2  => ipd.w[2] ,
          :s1  => ipd.σ[1] ,
          :s2  => ipd.σ[2] )
    else
        return OrderedDict(
          :w1  => ipd.w[1] ,
          :w2  => ipd.w[2] ,
          :w3  => ipd.w[3] ,
          :s1  => ipd.σ[1] ,
          :s2  => ipd.σ[2] ,
          :s3  => ipd.σ[3] )
    end
end
