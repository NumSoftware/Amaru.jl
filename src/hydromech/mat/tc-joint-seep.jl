 #This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export TCJointSeep

mutable struct TCJointSeepState<:IpState
    env::ModelEnv
    σ  ::Array{Float64,1} # stress
    w  ::Array{Float64,1} # relative displacements
    up::Float64          # effective plastic relative displacement
    Δλ ::Float64          # plastic multiplier
    h  ::Float64          # characteristic length from bulk elements
    uw ::Array{Float64,1} # interface pore pressure
    Vt ::Array{Float64,1} # transverse fluid velocity
    L  ::Array{Float64,1} 
    function TCJointSeepState(env::ModelEnv=ModelEnv())
        this = new(env)
        ndim = env.ndim
        this.σ   = zeros(ndim)
        this.w   = zeros(ndim)
        this.up  = 0.0
        this.Δλ  = 0.0
        this.h   = 0.0
        this.Vt  = zeros(2) 
        this.L   = zeros(ndim-1)
        this.uw  = zeros(3) # three layers
        return this
    end
end

mutable struct TCJointSeep<:Material
    E ::Float64       # Young's modulus
    ν ::Float64       # Poisson ratio
    ft::Float64       # tensile strength
    fc::Float64       # compressive strength
    ζ ::Float64       # factor ζ controls the elastic relative displacements (formerly α)
    wc::Float64       # critical crack opening
    ws::Float64       # openning at inflection (where the curve slope changes)
    softcurve::String # softening curve model ("linear" or bilinear" or "hordijk")
    α::Float64        # curvature coefficient
    γ::Float64        # factor for βres
    θ::Float64        # fator for the surface reduction speed
    βini::Float64     # initial curvature size
    γw::Float64       # specific weight of the fluid
    β ::Float64       # compressibility of fluid
    η ::Float64       # viscosity
    kt::Float64       # transverse leak-off coefficient
    w ::Float64       # initial fracture opening (longitudinal flow)
    fracture::Bool    # pre-existing fracture (true or false)

    function TCJointSeep(prms::Dict{Symbol,Float64})
        return  TCJointSeep(;prms...)
    end

     function TCJointSeep(;E=NaN, nu=NaN, ft=NaN, fc=NaN, zeta=NaN, alpha=1.0, gamma=1.0, theta=1.0, wc=NaN, ws=NaN, GF=NaN, Gf=NaN, softcurve="hordijk", gammaw=NaN, beta=0.0, eta=NaN, kt=NaN, w=0.0, fracture=false)  
        @check GF>0 || Gf>0
        @check E>0.0     
        @check 0<=nu<0.5 
        @check fc<0
        @check ft>=0     
        @check zeta>0    
        @check softcurve in ("linear", "bilinear", "hordijk")
        @check gammaw>0
        @check beta>= 0
        @check eta>=0  
        @check kt>=0   
        @check w>=0 
        
        if isnan(wc)
            if softcurve == "linear"
                wc = round(2*GF/ft, digits=10)
            elseif softcurve == "bilinear"
                if isnan(Gf)
                    wc = round(5*GF/ft, digits=10)
                    ws = round(wc*0.15, digits=10)
                else
                    wc = round((8*GF- 6*Gf)/ft, digits=10)
                    ws = round(1.5*Gf/ft, digits=10)
                end
            elseif softcurve == "hordijk"
                wc = round(GF/(0.1947019536*ft), digits=10)
            end
        end
        @check wc>0      
        @check isnan(ws) || ws>0
        
        a = (2*alpha*ft + alpha*fc - fc - √(alpha^2*fc^2 - 4*alpha^2*fc*ft + 4*alpha^2*ft^2 - 2*alpha*fc^2 + fc^2)) / (4*alpha-2)
        b = √(alpha*(2*a-fc)*(ft-a))
        βini = (b^2/ft^2)^alpha/(ft-a)

        this = new(E, nu, ft, fc, zeta, wc, ws, softcurve, alpha, gamma, theta, βini, gammaw, beta, eta, kt, w, fracture)
        return this
    end
end


# Type of corresponding state structure
compat_state_type(::Type{TCJointSeep}) = TCJointSeepState

# Element types that work with this material
compat_elem_types(::Type{TCJointSeep}) = (HMJoint,)


function mountD(mat::TCJointSeep, state::TCJointSeepState)
    if mat.fracture 
        # state.up = max(mat.wc, state.up)
        # state.w[1] = mat.w
        # mat.fracture = false
    end 

    invoke(mountD, Tuple{Material, IpState}, mat, state)
end


function update_state!(mat::TCJointSeep, state::TCJointSeepState, Δw::Array{Float64,1}, Δuw::Array{Float64,1},  G::Array{Float64,1}, BfUw::Array{Float64,1}, Δt::Float64)
    
    Δσ, status = invoke(update_state, Tuple{Material, IpState, Vector{Float64}}, mat, state, Δw)

    state.uw += Δuw
    state.Vt  = -mat.kt*G

    # crack opening
    if state.up == 0.0 || state.w[1] <= 0.0 
        w = 0.0
    else
        w = state.w[1]
    end

    state.L  =  ((w^3)/(12*mat.η))*BfUw

    return Δσ, state.Vt, state.L, status
end


function ip_state_vals(mat::TCJointSeep, state::TCJointSeepState)
    ndim = state.env.ndim
    if ndim == 3
       return OrderedDict(
          :jw1 => state.w[1],
          :jw2 => state.w[2],
          :jw3 => state.w[3],
          :js1 => state.σ[1],
          :js2 => state.σ[2],
          :js3 => state.σ[3],
          :jup => state.up,
          :juw => state.uw[3], # middle layer
          :vb  => state.Vt[1],
          :vt  => state.Vt[2])
    else
        return OrderedDict(
          :jw1 => state.w[1],
          :jw2 => state.w[2],
          :js1 => state.σ[1],
          :js2 => state.σ[2],
          :jup => state.up,
          :juw => state.uw[3], # middle layer
          :vb  => state.Vt[1],
          :vt  => state.Vt[2])
    end
end


function output_keys(mat::TCJointSeep)
    return Symbol[:jw1, :js1, :jup, :juw]
end
