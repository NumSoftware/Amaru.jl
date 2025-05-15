 #This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export MCJointSeep

mutable struct MCJointSeepState<:IpState
    ctx::Context
    σ  ::Array{Float64,1}  # stress
    w  ::Array{Float64,1}  # relative displacements
    Vt ::Array{Float64,1}  # transverse fluid velocity
    #D   ::Array{Float64,1}  # distance traveled by the fluid
    L  ::Array{Float64,1} 
    uw ::Array{Float64,1}  # interface pore pressure
    up ::Float64           # effective plastic relative displacement
    Δλ ::Float64           # plastic multiplier
    h  ::Float64           # characteristic length from bulk elements
    function MCJointSeepState(ctx::Context)
        this = new(ctx)
        ndim = ctx.ndim
        this.σ   = zeros(ndim)
        this.w   = zeros(ndim)
        this.Vt  = zeros(2) 
        this.L   = zeros(ndim-1)
        this.uw  = zeros(3) 
        this.up = 0.0
        this.Δλ  = 0.0
        this.h   = 0.0
        return this
    end
end


MCJointSeep_params = [
    FunInfo(:MCJointSeep, "Consitutive model for jointed rock with seepage"),
    KwArgInfo(:E, "Young's modulus", cond=:(E>0)),
    KwArgInfo(:nu, "Poisson ratio", cond=:(0<=nu<0.5)),
    KwArgInfo(:ft, "Tensile strength", cond=:(ft>0)),
    KwArgInfo(:mu, "Tangent of friction angle", cond=:(mu>=0)),
    KwArgInfo(:zeta, "Factor to control elastic relative displacements", cond=:(zeta>=0)),
    KwArgInfo(:wc, "Critical crack opening", 0.0, cond=:(wc>=0)),
    KwArgInfo(:GF, "Fracture energy", 0.0, cond=:(GF>=0)),
    KwArgInfo(:softmodel, "Softening model", :hordijk, values=(:linear, :bilinear, :hordijk, :soft, :custom), type=Symbol),
    KwArgInfo(:beta, "Compressibility of fluid", 0.0, cond=:(beta>=0)),
    KwArgInfo(:eta, "Viscosity", cond=:(eta>=0)),
    KwArgInfo(:kt, "Transverse leak-off coefficient", cond=:(kt>=0)),
    KwArgInfo(:w, "Initial fracture opening (longitudinal flow)", 0.0, cond=:(w>=0)),
    KwArgInfo(:fracture, "Pre-existing fracture", false, type=Bool),
]
@doc docstring(MCJointSeep_params) MCJointSeep

mutable struct MCJointSeep<:Material
    E  ::Float64       # Young's modulus
    ν  ::Float64       # Poisson ratio
    ft ::Float64     # tensile strength (internal variable)
    μ  ::Float64       # tangent of friction angle
    ζ  ::Float64       # factor ζ controls the elastic relative displacements (formerly α)
    wc ::Float64       # critical crack opening
    softmodel::Symbol  # softening curve model ("linear" or bilinear" or "hordijk")
    β  ::Float64       # compressibility of fluid
    η  ::Float64       # viscosity
    kt ::Float64       # transverse leak-off coefficient
    w ::Float64        # initial fracture opening (longitudinal flow)
    fracture::Bool     # pre-existing fracture (true or false)

     function MCJointSeep(; kwargs...)
        args = checkargs(kwargs, MCJointSeep_params)

        wc = args.wc
        softmodel = args.softmodel

        if wc==0.0
            GF = args.GF
            ft = args.ft
            if softmodel == :linear
                wc = round(2*GF/ft, sigdigits=5)
            elseif softmodel == :bilinear
                wc = round(5*GF/ft, sigdigits=5)
            elseif softmodel==:hordijk
                wc = round(GF/(0.1947019536*ft), sigdigits=5)  
            elseif softmodel==:soft
                wc = round(GF/(0.1947019536*ft), sigdigits=5)
            end
        end

        this = new(args.E, args.nu, args.ft, args.mu, args.zeta, wc, softmodel, args.beta, args.eta, args.kt, args.w, args.fracture)
        return this
    end
end


# Type of corresponding state structure
compat_state_type(::Type{MCJointSeep}, ::Type{HMJoint}, ctx::Context) = MCJointSeepState


function yield_func(mat::MCJointSeep, state::MCJointSeepState, σ::Array{Float64,1})
    ndim = state.ctx.ndim
    σmax = calc_σmax(mat, state, state.up)
    if ndim == 3
        return sqrt(σ[2]^2 + σ[3]^2) + (σ[1]-σmax)*mat.μ
    else
        return abs(σ[2]) + (σ[1]-σmax)*mat.μ
    end
end


function yield_deriv(mat::MCJointSeep, state::MCJointSeepState)
    ndim = state.ctx.ndim
    if ndim == 3
        return [ mat.μ, state.σ[2]/sqrt(state.σ[2]^2 + state.σ[3]^2), state.σ[3]/sqrt(state.σ[2]^2 + state.σ[3]^2)]
    else
        return [ mat.μ, sign(state.σ[2]) ]
    end
end


function potential_derivs(mat::MCJointSeep, state::MCJointSeepState, σ::Array{Float64,1})
    ndim = state.ctx.ndim
    if ndim == 3
            if σ[1] >= 0.0 
                # G1:
                r = [ 2.0*σ[1]*mat.μ^2, 2.0*σ[2], 2.0*σ[3]]
            else
                # G2:
                r = [ 0.0, 2.0*σ[2], 2.0*σ[3] ]
            end
    else
            if σ[1] >= 0.0 
                # G1:
                r = [ 2*σ[1]*mat.μ^2, 2*σ[2]]
            else
                # G2:
                r = [ 0.0, 2*σ[2] ]
            end
    end
    return r
end


function calc_σmax(mat::MCJointSeep, state::MCJointSeepState, up::Float64)
    if mat.softmodel == :linear
        if up < mat.wc
            a = mat.ft
            b = mat.ft/mat.wc
        else
            a = 0.0
            b = 0.0
        end
        σmax = a - b*up
    elseif mat.softmodel == :bilinear
        σs = 0.25*mat.ft
        if up < mat.ws
            a  = mat.ft 
            b  = (mat.ft - σs)/mat.ws
        elseif up < mat.wc
            a  = mat.wc*σs/(mat.wc-mat.ws)
            b  = σs/(mat.wc-mat.ws)
        else
            a = 0.0
            b = 0.0
        end
        σmax = a - b*up
    elseif mat.softmodel == :hordijk
        if up < mat.wc
            e = exp(1.0)
            z = (1 + 27*(up/mat.wc)^3)*e^(-6.93*up/mat.wc) - 28*(up/mat.wc)*e^(-6.93)
        else
            z = 0.0
        end
        σmax = z*mat.ft
    end
    return σmax
end


function σmax_deriv(mat::MCJointSeep, state::MCJointSeepState, up::Float64)
   # ∂σmax/∂up = dσmax
    if mat.softmodel == :linear
        if up < mat.wc
            b = mat.ft/mat.wc
        else
            b = 0.0
        end
        dσmax = -b
    elseif mat.softmodel == :bilinear
        σs = 0.25*mat.ft
        if up < mat.ws
            b  = (mat.ft - σs)/mat.ws
        elseif up < mat.wc
            b  = σs/(mat.wc-mat.ws)
        else
            b = 0.0
        end
        dσmax = -b
    elseif mat.softmodel == :hordijk
        if up < mat.wc
            e = exp(1.0)
            dz = ((81*up^2*e^(-6.93*up/mat.wc)/mat.wc^3) - (6.93*(1 + 27*up^3/mat.wc^3)*e^(-6.93*up/mat.wc)/mat.wc) - 0.02738402432/mat.wc)
        else
            dz = 0.0
        end
        dσmax = dz*mat.ft
    end
    return dσmax
end


function calc_kn_ks_De(mat::MCJointSeep, state::MCJointSeepState)
    ndim = state.ctx.ndim
    kn = mat.E*mat.ζ/state.h
    G  = mat.E/(2.0*(1.0+mat.ν))
    ks = G*mat.ζ/state.h

    if ndim == 3
        De = [  kn  0.0  0.0
               0.0   ks  0.0
               0.0  0.0   ks ]
    else
        De = [  kn   0.0
                 0.0  ks  ]
    end

    return kn, ks, De
end


function calc_Δλ(mat::MCJointSeep, state::MCJointSeepState, σtr::Array{Float64,1})
    ndim = state.ctx.ndim
    maxits = 100
    Δλ     = 0.0
    f      = 0.0
    up    = 0.0
    tol    = 1e-4      

    for i in 1:maxits
        μ      = mat.μ
        kn, ks, De = calc_kn_ks_De(mat, state)

        # quantities at n+1
        if ndim == 3
            if σtr[1]>0
                 σ     = [ σtr[1]/(1+2*Δλ*kn*μ^2),  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
                 dσdΔλ = [ -2*kn*μ^2*σtr[1]/(1+2*Δλ*kn*μ^2)^2,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2,  -2*ks*σtr[3]/(1+2*Δλ*ks)^2 ]
                 drdΔλ = [ -4*kn*μ^4*σtr[1]/(1+2*Δλ*kn*μ^2)^2,  -4*ks*σtr[2]/(1+2*Δλ*ks)^2,  -4*ks*σtr[3]/(1+2*Δλ*ks)^2 ]
            else
                 σ     = [ σtr[1],  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
                 dσdΔλ = [ 0,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2,  -2*ks*σtr[3]/(1+2*Δλ*ks)^2 ]
                 drdΔλ = [ 0,  -4*ks*σtr[2]/(1+2*Δλ*ks)^2,  -4*ks*σtr[3]/(1+2*Δλ*ks)^2 ]
            end
        else
            if σtr[1]>0
                 σ     = [ σtr[1]/(1+2*Δλ*kn*μ^2),  σtr[2]/(1+2*Δλ*ks) ]
                 dσdΔλ = [ -2*kn*μ^2*σtr[1]/(1+2*Δλ*kn*μ^2)^2,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2 ]
                 drdΔλ = [ -4*kn*μ^4*σtr[1]/(1+2*Δλ*kn*μ^2)^2,  -4*ks*σtr[2]/(1+2*Δλ*ks)^2 ]
            else
                 σ     = [ σtr[1],  σtr[2]/(1+2*Δλ*ks) ]
                 dσdΔλ = [ 0,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2 ]
                 drdΔλ = [ 0,  -4*ks*σtr[2]/(1+2*Δλ*ks)^2 ]
             end
        end
                 
         r      = potential_derivs(mat, state, σ)
         norm_r = norm(r)
         up    = state.up + Δλ*norm_r
         σmax   = calc_σmax(mat, state, up)
         m      = σmax_deriv(mat, state, up)
         dσmaxdΔλ = m*(norm_r + Δλ*dot(r/norm_r, drdΔλ))

        if ndim == 3
            f = sqrt(σ[2]^2 + σ[3]^2) + (σ[1]-σmax)*μ
            if (σ[2]==0 && σ[3]==0) 
                dfdΔλ = (dσdΔλ[1] - dσmaxdΔλ)*μ              
            else
                dfdΔλ = 1/sqrt(σ[2]^2 + σ[3]^2) * (σ[2]*dσdΔλ[2] + σ[3]*dσdΔλ[3]) + (dσdΔλ[1] - dσmaxdΔλ)*μ
            end
        else
            f = abs(σ[2]) + (σ[1]-σmax)*mat.μ
            dfdΔλ = sign(σ[2])*dσdΔλ[2] + (dσdΔλ[1] - dσmaxdΔλ)*μ
        end
        
        Δλ = Δλ - f/dfdΔλ

        abs(f) < tol && break

        if i == maxits || isnan(Δλ)
            warn("""PJoint: Could not find Δλ. This may happen when the system
            becomes hypostatic and thus the global stiffness matrix is nearly singular.
            Increasing the mesh refinement may result in a nonsingular matrix.
            """)
            warn("iterations=$i Δλ=$Δλ")
            return 0.0, failure("MCJointSeep: Could nof find Δλ.")
        end
    end
    return Δλ, success()
end


function calc_σ_upa(mat::MCJointSeep, state::MCJointSeepState, σtr::Array{Float64,1})
    ndim = state.ctx.ndim
    μ = mat.μ
    kn, ks, De = calc_kn_ks_De(mat, state)

    if ndim == 3
        if σtr[1] > 0
            σ = [σtr[1]/(1 + 2*state.Δλ*kn*(μ^2)), σtr[2]/(1 + 2*state.Δλ*ks), σtr[3]/(1 + 2*state.Δλ*ks)]
        else
            σ = [σtr[1], σtr[2]/(1 + 2*state.Δλ*ks), σtr[3]/(1 + 2*state.Δλ*ks)]
        end    
    else
        if σtr[1] > 0
            σ = [σtr[1]/(1 + 2*state.Δλ*kn*(μ^2)), σtr[2]/(1 + 2*state.Δλ*ks)]
        else
            σ = [σtr[1], σtr[2]/(1 + 2*state.Δλ*ks)]
        end    
    end
    state.σ = σ
    r = potential_derivs(mat, state, state.σ)
    state.up += state.Δλ*norm(r)
    return state.σ, state.up
end


function calcD(mat::MCJointSeep, state::MCJointSeepState)

    ndim = state.ctx.ndim
    kn, ks, De = calc_kn_ks_De(mat, state)

    if mat.fracture 
        state.up = mat.wc
    end 

    σmax = calc_σmax(mat, state, state.up)

    if state.Δλ == 0.0  # Elastic 
        return De
    elseif σmax == 0.0 
        Dep  = De*1e-10 
        return Dep
    else
        v    = yield_deriv(mat, state)
        r    = potential_derivs(mat, state, state.σ)
        y    = -mat.μ # ∂F/∂σmax
        m    = σmax_deriv(mat, state, state.up)  # ∂σmax/∂up

        #Dep  = De - De*r*v'*De/(v'*De*r - y*m*norm(r))

        if ndim == 3
            den = kn*r[1]*v[1] + ks*r[2]*v[2] + ks*r[3]*v[3] - y*m*norm(r)

            Dep = [   kn - kn^2*r[1]*v[1]/den    -kn*ks*r[1]*v[2]/den      -kn*ks*r[1]*v[3]/den
                     -kn*ks*r[2]*v[1]/den         ks - ks^2*r[2]*v[2]/den  -ks^2*r[2]*v[3]/den
                     -kn*ks*r[3]*v[1]/den        -ks^2*r[3]*v[2]/den        ks - ks^2*r[3]*v[3]/den ]
        else
            den = kn*r[1]*v[1] + ks*r[2]*v[2] - y*m*norm(r)

            Dep = [   kn - kn^2*r[1]*v[1]/den    -kn*ks*r[1]*v[2]/den      
                     -kn*ks*r[2]*v[1]/den         ks - ks^2*r[2]*v[2]/den  ]
        end

        return Dep
    end
end


function update_state!(mat::MCJointSeep, state::MCJointSeepState, Δw::Array{Float64,1}, Δuw::Array{Float64,1},  G::Array{Float64,1}, BfUw::Array{Float64,1}, Δt::Float64)
    ndim = state.ctx.ndim
    σini = copy(state.σ)

    kn, ks, De = calc_kn_ks_De(mat, state)

    if mat.fracture 
        state.up = mat.wc
    end 
    
    σmax = calc_σmax(mat, state, state.up) 

    if isnan(Δw[1]) || isnan(Δw[2])
        @warn "MCJointSeep: Invalid value for joint displacement: Δw = $Δw"
    end

    # σ trial and F trial
    σtr  = state.σ + De*Δw

    Ftr  = yield_func(mat, state, σtr) 

    # Elastic and EP integration
    if σmax == 0.0 && state.w[1] >= 0.0
        if ndim==3
            r1 = [ σtr[1]/kn, σtr[2]/ks, σtr[3]/ks ]
            r = r1/norm(r1)
            state.Δλ = norm(r1)
        else
            r1 = [ σtr[1]/kn, σtr[2]/ks ]
            r = r1/norm(r1)
            state.Δλ = norm(r1)  
        end

        state.up += state.Δλ
        state.σ = σtr - state.Δλ*De*r     

    elseif Ftr <= 0.0
        # Pure elastic increment
        state.Δλ = 0.0
        state.σ  = copy(σtr) 

    else    
        state.Δλ, status = calc_Δλ(mat, state, σtr)
        failed(status) && return state.σ,state.Vt,state.L, status
        state.σ, state.up = calc_σ_upa(mat, state, σtr)
                      
        # Return to surface:
        F  = yield_func(mat, state, state.σ)   
        if F > 1e-3
            @warn "MCJointSeep: Yield function value outside tolerance: $F"
        end
    end

    state.w += Δw
    Δσ = state.σ - σini

    state.uw += Δuw
    state.Vt  = -mat.kt*G
    #state.D  +=  state.Vt*Δt

    # compute crack aperture
    #if mat.w == 0.0
        if state.up == 0.0 || state.w[1] <= 0.0 
            w = 0.0
        else
            w = state.w[1]
        end
    #else
    #    if mat.w >= state.w[1]
    #        w = mat.w
    #    else 
    #        w = state.w[1]
    #    end
    #end 

    state.L  =  ((w^3)/(12*mat.η))*BfUw
    #state.S +=  state.L*Δt

    return Δσ, state.Vt, state.L, success()
end


function ip_state_vals(mat::MCJointSeep, state::MCJointSeepState)
    ndim = state.ctx.ndim
    if ndim == 3
       return OrderedDict(
          :jw   => state.w[1] ,
          :w2   => state.w[2] ,
          :w3   => state.w[3] ,
          :jσn   => state.σ[1] ,
          :s2   => state.σ[2] ,
          :s3   => state.σ[3] ,
          :up  => state.up  ,
          :uwf  => state.uw[3],
          :vb   => state.Vt[1],
          :vt   => state.Vt[2])
    else
        return OrderedDict(
          :jw   => state.w[1] ,
          :w2   => state.w[2] ,
          :jσn   => state.σ[1] ,
          :s2   => state.σ[2] ,
          :up  => state.up  ,
          :uwf  => state.uw[3],
          :vb   => state.Vt[1],
          :vt   => state.Vt[2])
    end
end
