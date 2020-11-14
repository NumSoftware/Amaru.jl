# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export JointLinSeep

mutable struct JointLinSeepIpState<:IpState
    env  ::ModelEnv
    V    ::Array{Float64,1} # fluid velocity
    #D    ::Array{Float64,1} # distance traveled by the fluid
    L    ::Array{Float64,1}
    #S    ::Array{Float64,1}
    uw   ::Array{Float64,1} # interface pore pressure
    h    ::Float64          # characteristic length from bulk elements
    function JointLinSeepIpState(env::ModelEnv=ModelEnv())
        this     = new(env)
        ndim     = env.ndim
        this.V   = zeros(2)
        #this.D   = zeros(2)
        this.L   = zeros(ndim-1)
        #this.S   = zeros(ndim-1)
        this.uw  = zeros(3)
        this.h   = 0.0
        return this
    end
end

mutable struct JointLinSeep<:Material
    γw ::Float64        # specific weight of the fluid
    β  ::Float64        # compressibility of fluid
    η  ::Float64        # viscosity
    kt ::Float64        # leak-off coefficient
    w ::Float64        # initial fracture opening (longitudinal flow)

    function JointLinSeep(prms::Dict{Symbol,Float64})
        return  JointLinSeep(;prms...)
    end

    function JointLinSeep(;gammaw=NaN, beta=0.0, eta=NaN, kt=NaN, w=0.0)
        
        gammaw>0    || error("Invalid value for gammaw: $gammaw")
        beta>= 0    || error("Invalid value for beta: $beta")
        eta>=0      || error("Invalid value for eta: $eta")
        kt>=0       || error("Invalid value for kt: $kt")
        w>=0       || error("Invalid value for w: $w")

        this = new(gammaw, beta, eta, kt, w)
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
    #ipd.D  +=  ipd.V*Δt
    ipd.L   =  ((mat.w^3)/(12*mat.η))*BfUw
    #ipd.S  +=  ipd.L*Δt

    return ipd.V, ipd.L
end

function ip_state_vals(mat::JointLinSeep, ipd::JointLinSeepIpState)
    return OrderedDict(
          :uwf => ipd.uw[3] ,
          :vb  => ipd.V[1] ,
          :vt  => ipd.V[2] )
end