 #This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export MCJoint

mutable struct MCJointState<:IpState
    env::ModelEnv
    σ  ::Array{Float64,1} # stress
    w  ::Array{Float64,1} # relative displacements
    up::Float64          # effective plastic relative displacement
    Δλ ::Float64          # plastic multiplier
    h  ::Float64          # characteristic length from bulk elements
    function MCJointState(env::ModelEnv=ModelEnv())
        this = new(env)
        ndim = env.ndim
        this.σ = zeros(ndim)
        this.w = zeros(ndim)
        this.up = 0.0
        this.Δλ  = 0.0
        this.h  = 0.0
        return this
    end
end

mutable struct MCJoint<:MatParams
    E  ::Float64      # Young's modulus
    ν  ::Float64      # Poisson ratio
    ft ::Float64      # tensile strength (internal variable)
    μ  ::Float64      # tangent of friction angle
    ζ  ::Float64      # factor ζ controls the elastic relative displacements (formerly α)
    wc ::Float64      # critical crack opening
    ws ::Float64      # openning at inflection (where the curve slope changes)
    softcurve::String # softening curve model ("linear" or bilinear" or "hordijk")

    function MCJoint(prms::Dict{Symbol,Float64})
        return  MCJoint(;prms...)
    end

    function MCJoint(;E=NaN, nu=NaN, ft=NaN, mu=NaN, zeta=NaN, wc=NaN, ws=NaN, GF=NaN, Gf=NaN, softcurve="hordijk")
        @check GF>0 || Gf>0 || wc>0
        @check E>0.0    
        @check 0<=nu<0.5
        @check ft>=0  
        @check mu>=0  
        @check zeta>0
        @check softcurve in ("linear", "bilinear", "hordijk", "soft")

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
            elseif softcurve=="hordijk" || softcurve=="soft"
                wc = round(GF/(0.1947019536*ft), digits=10) 
            end    
        end

        @check isnan(ws) || ws>0

        this = new(E, nu, ft, mu, zeta, wc, ws, softcurve)
        return this
    end
end



# Type of corresponding state structure
ip_state_type(matparams::MCJoint) = MCJointState


function yield_func(matparams::MCJoint, state::MCJointState, σ::Array{Float64,1})
    ndim = state.env.ndim
    σmax = calc_σmax(matparams, state, state.up)
    if ndim == 3
        return sqrt(σ[2]^2 + σ[3]^2) + (σ[1]-σmax)*matparams.μ
    else
        return abs(σ[2]) + (σ[1]-σmax)*matparams.μ
    end
end


function yield_deriv(matparams::MCJoint, state::MCJointState)
    ndim = state.env.ndim
    if ndim == 3
        return [ matparams.μ, state.σ[2]/sqrt(state.σ[2]^2 + state.σ[3]^2), state.σ[3]/sqrt(state.σ[2]^2 + state.σ[3]^2)]
    else
        return [ matparams.μ, sign(state.σ[2]) ]
    end
end


function potential_derivs(matparams::MCJoint, state::MCJointState, σ::Array{Float64,1})
    ndim = state.env.ndim
    if ndim == 3
        if σ[1] >= 0.0 
            # G1:
            r = [ 2.0*σ[1]*matparams.μ^2, 2.0*σ[2], 2.0*σ[3]]
        else
            # G2:
            r = [ 0.0, 2.0*σ[2], 2.0*σ[3] ]
        end
    else
        if σ[1] >= 0.0 
            # G1:
            r = [ 2*σ[1]*matparams.μ^2, 2*σ[2]]
        else
            # G2:
            r = [ 0.0, 2*σ[2] ]
        end
    end
    return r
end


function calc_σmax(matparams::MCJoint, state::MCJointState, up::Float64)
    if matparams.softcurve == "linear"
        if up < matparams.wc
            a = matparams.ft 
            b = matparams.ft /matparams.wc
        else
            a = 0.0
            b = 0.0
        end
        σmax = a - b*up
    elseif matparams.softcurve == "bilinear"
        σs = 0.25*matparams.ft 
        if up < matparams.ws
            a  = matparams.ft  
            b  = (matparams.ft  - σs)/matparams.ws
        elseif up < matparams.wc
            a  = matparams.wc*σs/(matparams.wc-matparams.ws)
            b  = σs/(matparams.wc-matparams.ws)
        else
            a = 0.0
            b = 0.0
        end
        σmax = a - b*up
    elseif matparams.softcurve == "hordijk"
        if up < matparams.wc
            e = exp(1.0)
            z = (1 + 27*(up/matparams.wc)^3)*e^(-6.93*up/matparams.wc) - 28*(up/matparams.wc)*e^(-6.93)           
        else
            z = 0.0
        end
        σmax = z*matparams.ft 
    elseif matparams.softcurve == "soft"
        m = 0.55
        a = 1.30837
        if up == 0.0
            z = 1.0
        elseif 0.0 < up < matparams.wc
            x = up/matparams.wc
            z = 1.0 - a^(1.0 - 1.0/x^m)
        else
            z = 0.0
        end
        σmax = z*matparams.ft
    end

    return σmax
end


function σmax_deriv(matparams::MCJoint, state::MCJointState, up::Float64)
    # ∂σmax/∂up = dσmax
    if matparams.softcurve == "linear"
        if up < matparams.wc
            b = matparams.ft /matparams.wc
        else
            b = 0.0
        end
        dσmax = -b
    elseif matparams.softcurve == "bilinear"
        σs = 0.25*matparams.ft 
        if up < matparams.ws
            b  = (matparams.ft  - σs)/matparams.ws
        elseif up < matparams.wc
            b  = σs/(matparams.wc-matparams.ws)
        else
            b = 0.0
        end
        dσmax = -b
    elseif matparams.softcurve == "hordijk"
        if up < matparams.wc
            e = exp(1.0)
            dz = ((81*up^2*e^(-6.93*up/matparams.wc)/matparams.wc^3) - (6.93*(1 + 27*up^3/matparams.wc^3)*e^(-6.93*up/matparams.wc)/matparams.wc) - 0.02738402432/matparams.wc)
        else
            dz = 0.0
        end
        dσmax = dz*matparams.ft 
    elseif matparams.softcurve == "soft"
        m = 0.55
        a = 1.30837

        if up == 0.0
            dz = 0.0
        elseif up < matparams.wc
            x = up/matparams.wc
            dz =  -m*log(a)*a^(1-x^-m)*x^(-m-1)/matparams.wc
        else
            dz = 0.0
        end
        dσmax = dz*matparams.ft 
    end

    return dσmax
end


function calc_kn_ks_De(matparams::MCJoint, state::MCJointState)
    ndim = state.env.ndim
    kn = matparams.E*matparams.ζ/state.h
    G  = matparams.E/(2.0*(1.0+matparams.ν))
    ks = G*matparams.ζ/state.h

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


function calc_Δλ(matparams::MCJoint, state::MCJointState, σtr::Array{Float64,1})
    ndim = state.env.ndim
    maxits = 100
    Δλ     = 0.0
    f      = 0.0
    up    = 0.0
    tol    = 1e-4      

    for i in 1:maxits
        μ      = matparams.μ
        kn, ks, De = calc_kn_ks_De(matparams, state)

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
                 
        r        = potential_derivs(matparams, state, σ)
        norm_r   = norm(r)
        up       = state.up + Δλ*norm_r
        σmax     = calc_σmax(matparams, state, up)
        m        = σmax_deriv(matparams, state, up)
        dσmaxdΔλ = m*(norm_r + Δλ*dot(r/norm_r, drdΔλ))

        if ndim == 3
            f = sqrt(σ[2]^2 + σ[3]^2) + (σ[1]-σmax)*μ
            if (σ[2]==0 && σ[3]==0) 
                dfdΔλ = (dσdΔλ[1] - dσmaxdΔλ)*μ
            else
                dfdΔλ = 1/sqrt(σ[2]^2 + σ[3]^2) * (σ[2]*dσdΔλ[2] + σ[3]*dσdΔλ[3]) + (dσdΔλ[1] - dσmaxdΔλ)*μ
            end
        else
            f = abs(σ[2]) + (σ[1]-σmax)*matparams.μ
            dfdΔλ = sign(σ[2])*dσdΔλ[2] + (dσdΔλ[1] - dσmaxdΔλ)*μ
        end

        Δλ = Δλ - f/dfdΔλ

        abs(f) < tol && break

        if i == maxits || isnan(Δλ)
            warn("""MCJoint: Could not find Δλ. This may happen when the system
            becomes hypostatic and thus the global stiffness matrix is nearly singular.
            Increasing the mesh refinement may result in a nonsingular matrix.
            """)
            warn("iterations=$i Δλ=$Δλ")
            return 0.0, failure("MCJoint: Could nof find Δλ.")
        end
    end
    return Δλ, success()
end


function calc_σ_upa(matparams::MCJoint, state::MCJointState, σtr::Array{Float64,1})
    ndim = state.env.ndim
    μ = matparams.μ
    kn, ks, De = calc_kn_ks_De(matparams, state)

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
    r = potential_derivs(matparams, state, state.σ)
    state.up += state.Δλ*norm(r)
    return state.σ, state.up
end


function mountD(matparams::MCJoint, state::MCJointState)
    ndim = state.env.ndim
    kn, ks, De = calc_kn_ks_De(matparams, state)
    σmax = calc_σmax(matparams, state, state.up)

    if state.Δλ == 0.0  # Elastic 
        return De
    elseif σmax == 0.0 
        # Dep  = De*1e-10 
        # Dep  = De*1e-5
        # Dep  = De*1e-4
        Dep  = De*1e-3
        return Dep
    else
        v    = yield_deriv(matparams, state)
        r    = potential_derivs(matparams, state, state.σ)
        y    = -matparams.μ # ∂F/∂σmax
        m    = σmax_deriv(matparams, state, state.up)  # ∂σmax/∂up

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


function update_state(matparams::MCJoint, state::MCJointState, Δw::Array{Float64,1})
    ndim = state.env.ndim
    σini = copy(state.σ)

    kn, ks, De = calc_kn_ks_De(matparams, state)
    σmax = calc_σmax(matparams, state, state.up)  

    if isnan(Δw[1]) || isnan(Δw[2])
        alert("MCJoint: Invalid value for joint displacement: Δw = $Δw")
    end

    # σ trial and F trial
    σtr  = state.σ + De*Δw

    Ftr  = yield_func(matparams, state, σtr) 

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
        state.Δλ, status = calc_Δλ(matparams, state, σtr)
        failed(status) && return state.σ, status

        state.σ, state.up = calc_σ_upa(matparams, state, σtr)
                      
        # Return to surface:
        F  = yield_func(matparams, state, state.σ)   
        F > 1e-3 && alert("MCJoint: Yield function value ($F) outside tolerance")

    end
    state.w += Δw
    Δσ = state.σ - σini
    return Δσ, success()
end


function ip_state_vals(matparams::MCJoint, state::MCJointState)
    ndim = state.env.ndim
    if ndim == 3
       return Dict(
          :jw1  => state.w[1] ,
          :jw2  => state.w[2] ,
          :jw3  => state.w[3] ,
          :js1  => state.σ[1] ,
          :js2  => state.σ[2] ,
          :js3  => state.σ[3] ,
          :jup => state.up
          )
    else
        return Dict(
          :jw1  => state.w[1] ,
          :jw2  => state.w[2] ,
          :js1  => state.σ[1] ,
          :js2  => state.σ[2] ,
          :jup => state.up
          )
    end
end


function output_keys(matparams::MCJoint)
    return Symbol[:jw1, :js1, :jup]
end
