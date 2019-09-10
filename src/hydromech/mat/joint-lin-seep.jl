# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export JointLinSeep

mutable struct JointLinSeepIpState<:IpState
    env  ::ModelEnv
    V    ::Array{Float64,1} # fluid velocity
    D    ::Array{Float64,1} # distance traveled by the fluid
    L    ::Array{Float64,1} 
    S    ::Array{Float64,1} 
    uw   ::Array{Float64,1} # interface pore pressure
    h    ::Float64          # characteristic length from bulk elements
    function JointLinSeepIpState(env::ModelEnv=ModelEnv())
        this     = new(env)
        ndim     = env.ndim
        this.V   = zeros(2) 
        this.D   = zeros(2) 
        this.L   = zeros(ndim-1)
        this.S   = zeros(ndim-1)
        this.uw  = zeros(3) 
        this.h   = 0.0
        return this
    end
end

mutable struct JointLinSeep<:Material
    k  ::Float64        # specific permeability
    γw ::Float64        # specific weight of the fluid
    α  ::Float64        # Biot's coefficient
    S  ::Float64        # Storativity coefficient
    β  ::Float64        # compressibility of fluid
    η  ::Float64        # viscosity
    kt ::Float64        # leak-off coefficient
    kl ::Float64        # initial fracture opening (longitudinal flow)

    function JointLinSeep(prms::Dict{Symbol,Float64})
        return  JointLinSeep(;prms...)
    end

    function JointLinSeep(;k=NaN, kappa=NaN, gammaw=NaN, alpha=NaN, S=NaN, n=NaN, Ks=NaN, Kw=NaN, beta=0.0, eta=NaN, kt=NaN, kl=0.0)

        !(isnan(kappa) || kappa>0) && error("Invalid value for kappa: $kappa")

        if isnan(k) 
            k = (kappa*gammaw)/eta # specific permeability = (intrinsic permeability * fluid specific weight)/viscosity
        end

        if isnan(S) 
            S = (alpha - n)/Ks + n/Kw # S = (alpha - porosity)/(bulk module of the solid) + (porosity)/(bulk module of the fluid) 
        end
        
        k>0         || error("Invalid value for k: $k")
        gammaw>0    || error("Invalid value for gammaw: $gammaw")
        0<alpha<=1.0|| error("Invalid value for alpha: $alpha")
        S>=0.0      || error("Invalid value for S: $S")
        beta>=0     || error("Invalid value for beta: $beta")
        eta>=0      || error("Invalid value for eta: $eta")
        kt>=0       || error("Invalid value for kt: $kt")
        kl>=0       || error("Invalid value for kl: $kl")

        this = new(k, gammaw, alpha, S, beta, eta, kt, kl)
        return this
    end
end

# Returns the element type that works with this material model
matching_elem_type(::JointLinSeep) = HydroJoint

# Type of corresponding state structure
ip_state_type(mat::JointLinSeep) = JointLinSeepIpState

function update_state!(mat::JointLinSeep, ipd::JointLinSeepIpState, Δuw::Array{Float64,1}, G::Array{Float64,1}, BfUw::Array{Float64,1}, Δt::Float64)
    ipd.uw +=  Δuw
    ipd.V   = -mat.kt*G
    ipd.D  +=  ipd.V*Δt
    ipd.L   =  ((mat.kl^3)/(12*mat.η))*BfUw
    ipd.S  +=  ipd.L*Δt
    
    return ipd.V, ipd.L
end

function ip_state_vals(mat::JointLinSeep, ipd::JointLinSeepIpState)
    
end
