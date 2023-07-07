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

mutable struct ElasticJointSeep<:MatParams
    E  ::Float64        # Young's modulus
    nu ::Float64        # Poisson ration 
    ζ  ::Float64        # factor ζ controls the elastic relative displacements 
    β  ::Float64        # compressibility of fluid
    η  ::Float64        # viscosity
    kt ::Float64        # leak-off coefficient
    w  ::Float64        # initial fracture opening (longitudinal flow)

    function ElasticJointSeep(prms::Dict{Symbol,Float64})
        return  ElasticJoint(;prms...)
    end

    function ElasticJointSeep(;E=NaN, nu=NaN, zeta=NaN, beta=0.0, eta=NaN, kt=NaN, w=0.0)
        @check E>0.0
        @check 0<=nu<0.5
        @check zeta>=0
        @check beta>= 0
        @check eta>=0
        @check kt>=0
        @check w>=0

        this = new(E, nu, zeta, beta, eta, kt, w)
        return this
    end
end



# Type of corresponding state structure
ip_state_type(::HydromechJointElem, ::ElasticJointSeep) = JointSeepState


function mountD(matparams::ElasticJointSeep, state::JointSeepState)
    ndim = state.env.ndim
    G  = matparams.E/(1.0+matparams.nu)/2.0
    kn = matparams.E*matparams.ζ/state.h
    ks =     G*matparams.ζ/state.h
    if ndim==2
        return  [  kn  0.0 
                  0.0   ks ]
    else
        return  [  kn  0.0  0.0
                  0.0   ks  0.0
                  0.0  0.0   ks ]
    end
end


function update_state(matparams::ElasticJointSeep, state::JointSeepState, Δu::Array{Float64,1}, Δuw::Array{Float64,1}, G::Array{Float64,1}, BfUw::Array{Float64,1}, Δt::Float64)
    ndim = state.env.ndim
    D  = mountD(matparams, state)
    Δσ = D*Δu

    state.w[1:ndim] += Δu
    state.σ[1:ndim] += Δσ

    state.uw += Δuw
    state.Vt = -matparams.kt*G
    #state.D +=  state.Vt*Δt

    # compute crack aperture
    if matparams.w == 0.0
        w = 0.0
    else
        if matparams.w >= state.w[1]
            w = matparams.w
        else 
            w = state.w[1]
        end
    end 

    state.L  =  ((w^3)/(12*matparams.η))*BfUw
    #state.S +=  state.L*Δt

    return Δσ, state.Vt, state.L, success()
end


function ip_state_vals(matparams::ElasticJointSeep, state::JointSeepState)
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
