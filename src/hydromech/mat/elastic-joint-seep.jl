# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticJointSeep

mutable struct JointSeepIpState<:IpState
    env::ModelEnv
    σ   ::Array{Float64,1}
    w   ::Array{Float64,1}
    uw  ::Float64  # fracture pore pressure
    upa ::Float64  # effective plastic relative displacement
    h   ::Float64
    time::Float64  # the time when the fracture opened
    function JointSeepIpState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ = zeros(3)
        this.w = zeros(3)
        this.uw = 0.0 
        this.upa = 0.0
        this.h = 0.0
        this.time = 0.0
        return this
    end
end

mutable struct ElasticJointSeep<:Material
    E::Float64 # Young modulus from bulk material
    nu::Float64 # Poisson ration from bulk material
    αi::Float64 # elastic displacement scale factor
    k ::Float64 # specific permeability
    γw::Float64 # water specific weight
    α ::Float64 # Biot's coefficient
    S ::Float64 # Storativity coefficient
    β  ::Float64  # coefficient of compressibility
    μ  ::Float64  # viscosity
    permeability::String # joint permeability specific permeability ("permeable" or "impermeable")

    function ElasticJointSeep(prms::Dict{Symbol,Float64})
        return  ElasticJoint(;prms...)
    end

    function ElasticJointSeep(;E=NaN, nu=NaN, alphai=NaN, k=NaN, gammaw=NaN, alpha=NaN, S=NaN, beta=0.0, mu=NaN, permeability="permeable")
        
        E<=0.0       && error("Invalid value for E: $E")
        !(0<=nu<0.5) && error("Invalid value for nu: $nu")
        mu<=0        && error("Invalid value for mu: $mu")
        !(permeability=="permeable" || permeability=="impermeable") && error("Invalid permeability: permeability must to be permeable or impermeable")

        this = new(E, nu, alphai, k, gammaw, alpha, S, beta, mu, permeability)
        return this
    end
end

# Returns the element type that works with this material model
matching_elem_type(::ElasticJointSeep) = HydroMechJoint

# Create a new instance of Ip data
new_ip_state(mat::ElasticJointSeep, env::ModelEnv) = JointSeepIpState(env)

function mountD(mat::ElasticJointSeep, ipd::JointSeepIpState)
    ndim = ipd.env.ndim
    G  = mat.E/(1.0+mat.nu)/2.0
    kn = mat.E*mat.α/ipd.h
    ks =     G*mat.α/ipd.h
    if ndim==2
        return [  kn  0.0 
                 0.0   ks ]
    else
        return  [  kn  0.0  0.0
                  0.0   ks  0.0
                  0.0  0.0   ks ]
    end
end

function stress_update(mat::ElasticJointSeep, ipd::JointSeepIpState, Δu::Array{Float64,1}, Δuw::Float64, time::Float64)
    ndim = ipd.env.ndim
    D  = mountD(mat, ipd)
    Δσ = D*Δu

    ipd.w[1:ndim] += Δu
    ipd.σ[1:ndim] += Δσ

    ipd.uw += Δuw

    return Δσ
end

function ip_state_vals(mat::ElasticJointSeep, ipd::JointSeepIpState)
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
