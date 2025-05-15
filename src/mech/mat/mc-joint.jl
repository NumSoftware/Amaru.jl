 #This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export MCJoint

mutable struct MCJointState<:IpState
    ctx::Context
    σ  ::Array{Float64,1} # stress
    w  ::Array{Float64,1} # relative displacements
    up::Float64          # effective plastic relative displacement
    Δλ ::Float64          # plastic multiplier
    h  ::Float64          # characteristic length from bulk elements
    function MCJointState(ctx::Context)
        this = new(ctx)
        ndim = ctx.ndim
        this.σ = zeros(ndim)
        this.w = zeros(ndim)
        this.up = 0.0
        this.Δλ  = 0.0
        this.h  = 0.0
        return this
    end
end


MCJoint_params = [
    FunInfo(:MCJoint, "Consitutive model for cohesive elements using the MC criterion"),
    KwArgInfo(:ft, "Tensile strength", cond=:(ft>0)),
    KwArgInfo(:mu, "Friction coefficient", cond=:(mu>=0)),
    KwArgInfo(:kn, "Normal stiffness per area", 0.0, cond=:(kn>0)),
    KwArgInfo(:ks, "Shear stiffness per area", 0.0, cond=:(ks>=0)),
    KwArgInfo(:E, "Young's modulus (when used as cohesive element)", 0.0, cond=:(E>=0)),
    KwArgInfo(:nu, "Poisson ratio (when used as cohesive element)", 0.0, cond=:(0<=nu<0.5)),
    KwArgInfo(:zeta, "Factor to control elastic relative displacements (when used as cohesive element)", 0.0, cond=:(zeta>=0)),
    KwArgInfo(:wc, "Critical crack opening", 0.0, cond=:(wc>=0)),
    KwArgInfo(:GF, "Fracture energy", 0.0, cond=:(GF>=0)),
    KwArgInfo(:softmodel, "Softening model", :hordijk, values=(:linear, :bilinear, :hordijk, :soft, :custom), type=Symbol),
]
@doc docstring(MCJoint_params) MCJoint

mutable struct MCJoint<:Material
    E  ::Float64      # Young's modulus
    ν  ::Float64      # Poisson ratio
    ft ::Float64      # tensile strength (internal variable)
    μ  ::Float64      # tangent of friction angle
    ζ  ::Float64      # factor ζ controls the elastic relative displacements (formerly α)
    kn ::Float64      # normal stiffness per area (optional, for fixed kn)
    ks ::Float64      # shear stiffness per area (optional, for fixed ks)
    wc ::Float64      # critical crack opening
    softmodel::Symbol    # softening curve model (:linear" or bilinear" or "hordijk"

    function MCJoint(; kwargs...)
        args = checkargs(kwargs, MCJoint_params)

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

        ζ = args.zeta
        kn, ks = args.kn, args.ks

        this = new(args.E, args.nu, args.ft, args.mu, ζ, kn, ks, wc, softmodel)
        return this
    end
end


# Type of corresponding state structure
compat_state_type(::Type{MCJoint}, ::Type{MechJoint}, ctx::Context) = MCJointState


function yield_func(mat::MCJoint, state::MCJointState, σ::Array{Float64,1})
    ndim = state.ctx.ndim
    σmax = calc_σmax(mat, state, state.up)
    if ndim == 3
        return sqrt(σ[2]^2 + σ[3]^2) + (σ[1]-σmax)*mat.μ
    else
        return abs(σ[2]) + (σ[1]-σmax)*mat.μ
    end
end


function yield_deriv(mat::MCJoint, state::MCJointState)
    ndim = state.ctx.ndim
    if ndim == 3
        return [ mat.μ, state.σ[2]/sqrt(state.σ[2]^2 + state.σ[3]^2), state.σ[3]/sqrt(state.σ[2]^2 + state.σ[3]^2)]
    else
        return [ mat.μ, sign(state.σ[2]) ]
    end
end


function potential_derivs(mat::MCJoint, state::MCJointState, σ::Array{Float64,1})
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


function calc_σmax(mat::MCJoint, state::MCJointState, up::Float64)
    if mat.softmodel == :linear
        if up < mat.wc
            a = mat.ft 
            b = mat.ft /mat.wc
        else
            a = 0.0
            b = 0.0
        end
        σmax = a - b*up
    elseif mat.softmodel == :bilinear
        σs = 0.25*mat.ft 
        if up < mat.ws
            a  = mat.ft  
            b  = (mat.ft  - σs)/mat.ws
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
    elseif mat.softmodel == :soft
        m = 0.55
        a = 1.30837
        if up == 0.0
            z = 1.0
        elseif 0.0 < up < mat.wc
            x = up/mat.wc
            z = 1.0 - a^(1.0 - 1.0/x^m)
        else
            z = 0.0
        end
        σmax = z*mat.ft
    end

    return σmax
end


function σmax_deriv(mat::MCJoint, state::MCJointState, up::Float64)
    # ∂σmax/∂up = dσmax
    if mat.softmodel == :linear
        if up < mat.wc
            b = mat.ft /mat.wc
        else
            b = 0.0
        end
        dσmax = -b
    elseif mat.softmodel == :bilinear
        σs = 0.25*mat.ft 
        if up < mat.ws
            b  = (mat.ft  - σs)/mat.ws
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
    elseif mat.softmodel == :soft
        m = 0.55
        a = 1.30837

        if up == 0.0
            dz = 0.0
        elseif up < mat.wc
            x = up/mat.wc
            dz =  -m*log(a)*a^(1-x^-m)*x^(-m-1)/mat.wc
        else
            dz = 0.0
        end
        dσmax = dz*mat.ft 
    end

    return dσmax
end


function calc_kn_ks_De(mat::MCJoint, state::MCJointState)
    if mat.kn>0 && mat.ks>0
        kn = mat.kn
        ks = mat.ks
    else
        kn = mat.E*mat.ζ/state.h
        G  = mat.E/(2*(1 + mat.ν))
        ks = G*mat.ζ/state.h
    end
    
    ndim = state.ctx.ndim
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


function calc_Δλ(mat::MCJoint, state::MCJointState, σtr::Array{Float64,1})
    ndim   = state.ctx.ndim
    maxits = 100
    Δλ     = 0.0
    f      = 0.0
    up     = 0.0
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
                 
        r        = potential_derivs(mat, state, σ)
        norm_r   = norm(r)
        up       = state.up + Δλ*norm_r
        σmax     = calc_σmax(mat, state, up)
        m        = σmax_deriv(mat, state, up)
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
            # warn("""MCJoint: Could not find Δλ. This may happen when the system
            # becomes hypostatic and thus the global stiffness matrix is nearly singular.
            # Increasing the mesh refinement may result in a nonsingular matrix.
            # """)
            # warn("iterations=$i Δλ=$Δλ")
            return 0.0, failure("MCJoint: Could nof find Δλ.")
        end
    end
    return Δλ, success()
end


function calc_σ_upa(mat::MCJoint, state::MCJointState, σtr::Array{Float64,1})
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


function calcD(mat::MCJoint, state::MCJointState)
    ndim = state.ctx.ndim
    kn, ks, De = calc_kn_ks_De(mat, state)
    σmax = calc_σmax(mat, state, state.up)

    if state.Δλ == 0.0  # Elastic 
        return De
    elseif σmax == 0.0 
        # Dep  = De*1e-10 
        # Dep  = De*1e-5
        # Dep  = De*1e-4
        Dep  = De*1e-3
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


function update_state!(mat::MCJoint, state::MCJointState, Δw::Array{Float64,1})
    ndim = state.ctx.ndim
    σini = copy(state.σ)

    kn, ks, De = calc_kn_ks_De(mat, state)
    σmax = calc_σmax(mat, state, state.up)  

    if isnan(Δw[1]) || isnan(Δw[2])
        alert("MCJoint: Invalid value for joint displacement: Δw = $Δw")
    end

    # σ trial and F trial
    σtr  = state.σ + De*Δw

    Ftr  = yield_func(mat, state, σtr) 

    # Elastic and EP integration
    if σmax == 0.0 && state.w[1] >= 0.0
        # Return to apex:
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
        # Plastic increment
        state.Δλ, status = calc_Δλ(mat, state, σtr)
        failed(status) && return state.σ, status

        state.σ, state.up = calc_σ_upa(mat, state, σtr)
                      
        # Return to surface:
        F  = yield_func(mat, state, state.σ)   
        F > 1e-3 && alert("MCJoint: Yield function value ($F) outside tolerance")

    end
    state.w += Δw
    Δσ = state.σ - σini
    return Δσ, success()
end


function ip_state_vals(mat::MCJoint, state::MCJointState)
    ndim = state.ctx.ndim
    if ndim == 3
       return Dict(
          :jw  => state.w[1] ,
          :jw2  => state.w[2] ,
          :jw3  => state.w[3] ,
          :jσn  => state.σ[1] ,
          :js2  => state.σ[2] ,
          :js3  => state.σ[3] ,
          :jup => state.up
          )
    else
        return Dict(
          :jw  => state.w[1] ,
          :jw2  => state.w[2] ,
          :jσn  => state.σ[1] ,
          :js2  => state.σ[2] ,
          :jup => state.up
          )
    end
end


function output_keys(mat::MCJoint)
    return Symbol[:jw, :jσn, :jup]
end
