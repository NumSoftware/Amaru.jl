# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticJointSeep

mutable struct JointSeepState<:IpState
    env  ::ModelEnv
    σ    ::Array{Float64,1} # stress
    w    ::Array{Float64,1} # relative displacements
    Vt   ::Array{Float64,1} # transverse fluid velocity
    #D    ::Array{Float64,1} # distance traveled by the fluid
    L    ::Array{Float64,1} 
    #S    ::Array{Float64,1}
    uw   ::Array{Float64,1} # interface pore pressure
    h    ::Float64          # characteristic length from bulk elements
    up  ::Float64          # effective plastic relative displacement
    function JointSeepState(env::ModelEnv=ModelEnv())
        this     = new(env)
        ndim     = env.ndim
        this.σ   = zeros(3)
        this.w   = zeros(3)
        this.Vt  = zeros(2) 
        #this.D   = zeros(2) 
        this.L   = zeros(ndim-1)
        #this.S   = zeros(ndim-1)
        this.uw  = zeros(3) 
        this.h   = 0.0
        this.up = 0.0
        return this
    end
end

mutable struct ElasticJointSeep<:Material
    E  ::Float64        # Young's modulus
    nu ::Float64        # Poisson ration 
    ζ  ::Float64        # factor ζ controls the elastic relative displacements 
    γw ::Float64        # specific weight of the fluid
    β  ::Float64        # compressibility of fluid
    η  ::Float64        # viscosity
    kt ::Float64        # leak-off coefficient
    w  ::Float64        # initial fracture opening (longitudinal flow)

    function ElasticJointSeep(prms::Dict{Symbol,Float64})
        return  ElasticJoint(;prms...)
    end

    function ElasticJointSeep(;E=NaN, nu=NaN, zeta=NaN, gammaw=NaN, beta=0.0, eta=NaN, kt=NaN, w=0.0)
    
        E>0.0       || error("Invalid value for E: $E")
        0<=nu<0.5   || error("Invalid value for nu: $nu") 
        zeta>=0     || error("Invalid value for zeta: $zeta")
        gammaw>0    || error("Invalid value for gammaw: $gammaw")
        beta>= 0    || error("Invalid value for beta: $beta")
        eta>=0      || error("Invalid value for eta: $eta")
        kt>=0       || error("Invalid value for kt: $kt")
        w>=0        || error("Invalid value for w: $w")

        this = new(E, nu, zeta, gammaw, beta, eta, kt, w)
        return this
    end
end

# Returns the element type that works with this material model
matching_elem_type(::ElasticJointSeep) = HydroMechJoint

# Type of corresponding state structure
ip_state_type(mat::ElasticJointSeep) = JointSeepState

function mountD(mat::ElasticJointSeep, state::JointSeepState)
    ndim = state.env.ndim
    G  = mat.E/(1.0+mat.nu)/2.0
    kn = mat.E*mat.ζ/state.h
    ks =     G*mat.ζ/state.h
    if ndim==2
        return  [  kn  0.0 
                  0.0   ks ]
    else
        return  [  kn  0.0  0.0
                  0.0   ks  0.0
                  0.0  0.0   ks ]
    end
end

function stress_update(mat::ElasticJointSeep, state::JointSeepState, Δu::Array{Float64,1}, Δuw::Array{Float64,1}, G::Array{Float64,1}, BfUw::Array{Float64,1}, Δt::Float64)
    ndim = state.env.ndim
    D  = mountD(mat, state)
    Δσ = D*Δu

    state.w[1:ndim] += Δu
    state.σ[1:ndim] += Δσ

    state.uw += Δuw
    state.Vt = -mat.kt*G
    #state.D +=  state.Vt*Δt

    # compute crack aperture
    if mat.w == 0.0
        w = 0.0
    else
        if mat.w >= state.w[1]
            w = mat.w
        else 
            w = state.w[1]
        end
    end 

    state.L  =  ((w^3)/(12*mat.η))*BfUw
    #state.S +=  state.L*Δt

    return Δσ, state.Vt, state.L, success()
end

function ip_state_vals(mat::ElasticJointSeep, state::JointSeepState)
    ndim = state.env.ndim
    if ndim == 2
        return OrderedDict(
          :w1  => state.w[1]  ,
          :w2  => state.w[2]  ,
          :s1  => state.σ[1]  ,
          :s2  => state.σ[2]  ,
          :uwf => state.uw[3] ,
          :vb  => state.Vt[1] ,
          :vt  => state.Vt[2] )
    else
        return OrderedDict(
          :w1  => state.w[1]  ,
          :w2  => state.w[2]  ,
          :w3  => state.w[3]  ,
          :s1  => state.σ[1]  ,
          :s2  => state.σ[2]  ,
          :s3  => state.σ[3]  ,
          :uwf => state.uw[3] ,
          :vb  => state.Vt[1] ,
          :vt  => state.Vt[2] )
    end
end
