# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticJointSeep

mutable struct JointSeepIpState<:IpState
    env  ::ModelEnv
    σ    ::Array{Float64,1} # stress
    w    ::Array{Float64,1} # relative displacements
    uw   ::Array{Float64,1} # interface pore pressure
    G    ::Array{Float64,1} # longitudinal flow gradient 
    h    ::Float64          # characteristic length from bulk elements
    t    ::Float64          # time when the fracture opened
    function JointSeepIpState(env::ModelEnv=ModelEnv())
        this = new(env)
        ndim = env.ndim
        this.σ = zeros(3)
        this.w = zeros(3)
        this.uw = zeros(3) 
        this.G  = zeros(ndim-1)
        this.h = 0.0
        this.t = 0.0
        return this
    end
end

mutable struct ElasticJointSeep<:Material
    E ::Float64        # Young's modulus
    nu::Float64        # Poisson ration 
    ζ ::Float64        # factor ζ controls the elastic relative displacements 
    k ::Float64        # specific permeability
    γw::Float64        # specific weight of the fluid
    α ::Float64        # Biot's coefficient
    S ::Float64        # Storativity coefficient
    β ::Float64        # compressibility of fluid
    η ::Float64        # viscosity
    kt::Float64        # leak-off coefficient
    kl ::Float64       # longitudinal permeability coefficient
    permeability::Bool # joint permeability ("true" or "false")

    function ElasticJointSeep(prms::Dict{Symbol,Float64})
        return  ElasticJoint(;prms...)
    end

    function ElasticJointSeep(;E=NaN, nu=NaN, zeta=NaN, k=NaN, kappa=NaN, gammaw=NaN, alpha=NaN, S=NaN, n=NaN, Ks=NaN, Kw=NaN, beta=0.0, eta=NaN, kt=NaN, kl=0.0, permeability=true)

        !(isnan(kappa) || kappa>0) && error("Invalid value for kappa: $kappa")

        if isnan(k) 
            k = (kappa*gammaw)/eta # specific permeability = (intrinsic permeability * fluid specific weight)/viscosity
        end

        if isnan(S) 
            S = (alpha - n)/Ks + n/Kw # S = (alpha - porosity)/(bulk module of the solid) + (porosity)/(bulk module of the fluid) 
        end
        
        E>0.0       || error("Invalid value for E: $E")
        0<=nu<0.5   || error("Invalid value for nu: $nu") 
        zeta>=0     || error("Invalid value for zeta: $zeta")
        k>0         || error("Invalid value for k: $k")
        gammaw>0    || error("Invalid value for gammaw: $gammaw")
        0<alpha<=1.0|| error("Invalid value for alpha: $alpha")
        S>=0.0      || error("Invalid value for S: $S")
        beta>=0     || error("Invalid value for beta: $beta")
        eta>=0      || error("Invalid value for eta: $eta")
        kt>=0       || error("Invalid value for kt: $kt")
        kl>=0       || error("Invalid value for kl: $kl")
        (permeability==true || permeability==false) || error("Invalid permeability: permeability must to be true or false")

        this = new(E, nu, zeta, k, gammaw, alpha, S, beta, eta, kt, kl, permeability)
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
    kn = mat.E*mat.ζ/ipd.h
    ks =     G*mat.ζ/ipd.h
    if ndim==2
        return  [  kn  0.0 
                  0.0   ks ]
    else
        return  [  kn  0.0  0.0
                  0.0   ks  0.0
                  0.0  0.0   ks ]
    end
end

function stress_update(mat::ElasticJointSeep, ipd::JointSeepIpState, Δu::Array{Float64,1}, Δuw::Array{Float64,1}, G::Array{Float64,1})
    ndim = ipd.env.ndim
    D  = mountD(mat, ipd)
    Δσ = D*Δu

    ipd.w[1:ndim] += Δu
    ipd.σ[1:ndim] += Δσ

    ipd.uw += Δuw
    ipd.G   = G

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
