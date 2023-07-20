# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ConstPermeabilityJoint

mutable struct ConstPermeabilityJointState<:IpState
    env  ::ModelEnv
    V    ::Array{Float64,1} # fluid velocity
    #D    ::Array{Float64,1} # distance traveled by the fluid
    L    ::Array{Float64,1}
    #S    ::Array{Float64,1}
    uw   ::Array{Float64,1} # interface pore pressure
    h    ::Float64          # characteristic length from bulk elements
    function ConstPermeabilityJointState(env::ModelEnv=ModelEnv())
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

mutable struct ConstPermeabilityJoint<:Material
    γw ::Float64        # specific weight of the fluid
    β  ::Float64        # compressibility of fluid
    η  ::Float64        # viscosity
    kt ::Float64        # leak-off coefficient
    w ::Float64        # initial fracture opening (longitudinal flow)

    function ConstPermeabilityJoint(;gammaw=NaN, beta=0.0, eta=NaN, kt=NaN, w=0.0)
        names = (E="Young modulus", nu="Poisson ratio")
        required = (:E, :nu)
        @checkmissing params required names

        params  = values(params)
        E       = params.E
        nu      = params.nu

        @check E>=0.0
        @check 0<=nu<0.5
        return new(E, nu)
        return this
        
        @check gammaw>0
        @check beta>=0
        @check eta>=0  
        @check kt>=0   
        @check w>=0    

        this = new(gammaw, beta, eta, kt, w)
        return this
    end
end


# Type of corresponding state structure
ip_state_type(::Type{ConstPermeabilityJoint}) = ConstPermeabilityJointState

# Element types that work with this material
matching_elem_types(::Type{ConstPermeabilityJoint}) = (HydroJoint,)


function update_state!(mat::ConstPermeabilityJoint, state::ConstPermeabilityJointState, Δuw::Array{Float64,1}, G::Array{Float64,1}, BfUw::Array{Float64,1}, Δt::Float64)
    state.uw +=  Δuw
    state.V   = -mat.kt*G
    #state.D  +=  state.V*Δt
    state.L   =  ((mat.w^3)/(12*mat.η))*BfUw
    #state.S  +=  state.L*Δt

    return state.V, state.L
end


function ip_state_vals(mat::ConstPermeabilityJoint, state::ConstPermeabilityJointState)
    return OrderedDict(
          :uwf => state.uw[3] ,
          :vb  => state.V[1] ,
          :vt  => state.V[2] )
end