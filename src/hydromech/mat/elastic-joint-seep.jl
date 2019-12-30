# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticJointSeep

mutable struct JointSeepIpState<:IpState
    env  ::ModelEnv
    σ    ::Array{Float64,1} # stress
    w    ::Array{Float64,1} # relative displacements
    Vt   ::Array{Float64,1} # transverse fluid velocity
    D    ::Array{Float64,1} # distance traveled by the fluid
    L    ::Array{Float64,1} 
    S    ::Array{Float64,1}
    uw   ::Array{Float64,1} # interface pore pressure
    h    ::Float64          # characteristic length from bulk elements
    upa  ::Float64          # effective plastic relative displacement
    function JointSeepIpState(env::ModelEnv=ModelEnv())
        this     = new(env)
        ndim     = env.ndim
        this.σ   = zeros(3)
        this.w   = zeros(3)
        this.Vt  = zeros(2) 
        this.D   = zeros(2) 
        this.L   = zeros(ndim-1)
        this.S   = zeros(ndim-1)
        this.uw  = zeros(3) 
        this.h   = 0.0
        this.upa = 0.0
        return this
    end
end

mutable struct ElasticJointSeep<:Material
    E  ::Float64        # Young's modulus
    nu ::Float64        # Poisson ration 
    ζ  ::Float64        # factor ζ controls the elastic relative displacements 
    k  ::Float64        # specific permeability
    γw ::Float64        # specific weight of the fluid
    α  ::Float64        # Biot's coefficient
    S  ::Float64        # Storativity coefficient
    β  ::Float64        # compressibility of fluid
    η  ::Float64        # viscosity
    kt ::Float64        # leak-off coefficient
    kl ::Float64        # initial fracture opening (longitudinal flow)

    function ElasticJointSeep(prms::Dict{Symbol,Float64})
        return  ElasticJoint(;prms...)
    end

    function ElasticJointSeep(;E=NaN, nu=NaN, zeta=NaN, k=NaN, kappa=NaN, gammaw=NaN, alpha=NaN, S=NaN, n=NaN, Ks=NaN, Kw=NaN, beta=0.0, eta=NaN, kt=NaN, kl=0.0)

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

        this = new(E, nu, zeta, k, gammaw, alpha, S, beta, eta, kt, kl)
        return this
    end
end

# Returns the element type that works with this material model
matching_elem_type(::ElasticJointSeep) = HydroMechJoint

# Type of corresponding state structure
ip_state_type(mat::ElasticJointSeep) = JointSeepIpState

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

function stress_update(mat::ElasticJointSeep, ipd::JointSeepIpState, Δu::Array{Float64,1}, Δuw::Array{Float64,1}, G::Array{Float64,1}, BfUw::Array{Float64,1}, Δt::Float64)
    ndim = ipd.env.ndim
    D  = mountD(mat, ipd)
    Δσ = D*Δu

    ipd.w[1:ndim] += Δu
    ipd.σ[1:ndim] += Δσ

    ipd.uw += Δuw
    ipd.Vt = -mat.kt*G
    ipd.D +=  ipd.Vt*Δt

    # compute crack aperture
    if mat.kl == 0.0
        kl = 0.0
    else
        if mat.kl >= ipd.w[1]
            kl = mat.kl
        else 
            kl = ipd.w[1]
        end
    end 

    ipd.L  =  ((kl^3)/(12*mat.η))*BfUw
    ipd.S +=  ipd.L*Δt

    return Δσ, ipd.Vt, ipd.L
end

function ip_state_vals(mat::ElasticJointSeep, ipd::JointSeepIpState)
    ndim = ipd.env.ndim
    if ndim == 2
        return OrderedDict(
          :w1  => ipd.w[1]  ,
          :w2  => ipd.w[2]  ,
          :s1  => ipd.σ[1]  ,
          :s2  => ipd.σ[2]  ,
          :uw  => ipd.uw    ,
          :vb  => ipd.Vt[1] ,
          :vt  => ipd.Vt[2] )
    else
        return OrderedDict(
          :w1  => ipd.w[1]  ,
          :w2  => ipd.w[2]  ,
          :w3  => ipd.w[3]  ,
          :s1  => ipd.σ[1]  ,
          :s2  => ipd.σ[2]  ,
          :s3  => ipd.σ[3]  ,
          :uw  => ipd.uw    ,
          :vb  => ipd.Vt[1] ,
          :vt  => ipd.Vt[2] )
    end
end
