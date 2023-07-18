 #This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export PJoint

mutable struct PJointState<:IpState
    env::ModelEnv
    σ  ::Array{Float64,1} # stress
    w  ::Array{Float64,1} # relative displacements
    up::Float64          # effective plastic relative displacement
    Δλ ::Float64          # plastic multiplier
    h  ::Float64          # characteristic length from bulk elements
    function PJointState(env::ModelEnv=ModelEnv())
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

mutable struct PJoint<:Material
    E ::Float64      # Young's modulus
    ν ::Float64      # Poisson ratio
    ft::Float64      # tensile strength (internal variable)
    fc::Float64      # compressive strength (internal variable)
    ζ ::Float64      # factor ζ controls the elastic relative displacements (formerly α)
    wc::Float64      # critical crack opening
    ws::Float64      # openning at inflection (where the curve slope changes)
    softcurve::String # softening curve model ("linear" or bilinear" or "hordijk")
    α::Float64 # fator for reducing the curve
    γ::Float64 # factor for βres

    function PJoint(prms::Dict{Symbol,Float64})
        return  PJoint(;prms...)
    end

    function PJoint(;E=NaN, nu=NaN, fc=NaN, ft=NaN, zeta=NaN, wc=NaN, ws=NaN, GF=NaN, Gf=NaN, softcurve="bilinear")

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

        E>0.0      || error("Invalid value for E: $E")
        0<=nu<0.5  || error("Invalid value for nu: $nu")
        @check fc<0  
        @check ft>=0  
        zeta>0     || error("Invalid value for zeta: $zeta")
        wc>0       || error("Invalid value for wc: $wc")
        (isnan(ws) || ws>0) || error("Invalid value for ws: $ws")
        softcurve in ("linear", "bilinear", "hordijk") || error("Invalid softcurve: softcurve must to be linear or bilinear or hordijk")

        this = new(E, nu, ft, fc, zeta, wc, ws, softcurve, 1.5, 0.1)
        return this
    end
end



# Type of corresponding state structure
ip_state_type(mat::PJoint) = PJointState


function yield_func(mat::PJoint, state::PJointState, σ::Array{Float64,1}, σmax::Float64)
    ndim = state.env.ndim
    fc, ft = mat.fc, mat.ft

    βini = 2*ft - fc -2*√(ft^2 - fc*ft)
    βres = mat.γ*βini
    β = βres + (βini-βres)*(σmax/ft)^mat.α
    if ndim == 3
        return β*(σ[1] - σmax) + σ[2]^2 + σ[3]^2
    else
        return β*(σ[1] - σmax) + σ[2]^2
    end
end


function yield_derivs(mat::PJoint, state::PJointState, σ::Array{Float64,1}, σmax::Float64)
    ndim = state.env.ndim
    fc, ft = mat.fc, mat.ft
    βini = 2*ft - fc -2*√(ft^2 - fc*ft)
    βres = mat.γ*βini
    β = βres + (βini-βres)*(σmax/ft)^mat.α

    if ndim == 3
        return [ β, 2*σ[2], 2*σ[3] ]
    else
        return [ β, 2*σ[2] ]
    end
end


function potential_derivs(mat::PJoint, state::PJointState, σ::Array{Float64,1})
    ndim = state.env.ndim
    if ndim == 3
        if σ[1] > 0.0 
            # G1:
            r = [ 2.0*σ[1], 2.0*σ[2], 2.0*σ[3]]
        else
            # G2:
            r = [ 0.0, 2.0*σ[2], 2.0*σ[3] ]
        end
    else
        if σ[1] > 0.0 
            # G1:
            r = [ 2*σ[1], 2*σ[2]]
        else
            # G2:
            r = [ 0.0, 2*σ[2] ]
        end
    end
    return r
end


function calc_σmax(mat::PJoint, state::PJointState, up::Float64)
    if mat.softcurve == "linear"
        if up < mat.wc
            a = mat.ft 
            b = mat.ft /mat.wc
        else
            a = 0.0
            b = 0.0
        end
        σmax = a - b*up
    elseif mat.softcurve == "bilinear"
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
    elseif mat.softcurve == "hordijk"
        if up < mat.wc
            e = exp(1.0)
            z = (1 + 27*(up/mat.wc)^3)*e^(-6.93*up/mat.wc) - 28*(up/mat.wc)*e^(-6.93)
        else
            z = 0.0
        end
        σmax = z*mat.ft 
    end

    # σmax<0.001*mat.ft  && (σmax=0.0)

    return σmax
end


function deriv_σmax_upa(mat::PJoint, state::PJointState, up::Float64)
    # ∂σmax/∂up = dσmax
    if mat.softcurve == "linear"
        if up < mat.wc
            b = mat.ft /mat.wc
        else
            b = 0.0
        end
        dσmax = -b
    elseif mat.softcurve == "bilinear"
        σs = 0.25*mat.ft 
        if up < mat.ws
            b  = (mat.ft  - σs)/mat.ws
        elseif up < mat.wc
            b  = σs/(mat.wc-mat.ws)
        else
            b = 0.0
        end
        dσmax = -b
    elseif mat.softcurve == "hordijk"
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


function calc_kn_ks(mat::PJoint, state::PJointState)
    kn = mat.E*mat.ζ/state.h
    G  = mat.E/(2.0*(1.0+mat.ν))
    ks = G*mat.ζ/state.h

    return kn, ks
end


function calc_Δλ(mat::PJoint, state::PJointState, σtr::Array{Float64,1})
    ndim = state.env.ndim
    maxits = 20
    Δλ     = 0.0
    f      = 0.0
    up    = 0.0
    tol    = 1e-2
    fc, ft = mat.fc, mat.ft
    βini = 2*ft - fc -2*√(ft^2 - fc*ft)
    βres = mat.γ*βini
    α    = mat.α
    
    nits = 0

    for i in 1:maxits
        nits = i
        kn, ks = calc_kn_ks(mat, state)

        # quantities at n+1
        if ndim == 3
            if σtr[1]>0
                 σ     = [ σtr[1]/(1+2*Δλ*kn),  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
                 dσdΔλ = [ -2*kn*σtr[1]/(1+2*Δλ*kn)^2,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2,  -2*ks*σtr[3]/(1+2*Δλ*ks)^2 ]
            else
                 σ     = [ σtr[1],  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
                 dσdΔλ = [ 0,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2,  -2*ks*σtr[3]/(1+2*Δλ*ks)^2 ]
            end
        else
            if σtr[1]>0
                 σ     = [ σtr[1]/(1+2*Δλ*kn),  σtr[2]/(1+2*Δλ*ks) ]
                 dσdΔλ = [ -2*kn*σtr[1]/(1+2*Δλ*kn)^2,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2 ]
            else
                 σ     = [ σtr[1],  σtr[2]/(1+2*Δλ*ks) ]
                 dσdΔλ = [ 0,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2 ]
             end
        end

        drdΔλ = 2*dσdΔλ
                 
        r      = potential_derivs(mat, state, σ)
        norm_r = norm(r)
        up    = state.up + Δλ*norm_r
        σmax   = calc_σmax(mat, state, up)
        # β      = βres + (βini-βres)/ft*σmax
        β      = βres + (βini-βres)*(σmax/ft)^mat.α

        if ndim == 3
            f = β*(σ[1] - σmax) + σ[2]^2 + σ[3]^2
            dfdσ = [ β, 2*σ[2], 2*σ[3] ]
        else
            f = β*(σ[1] - σmax) + σ[2]^2
            dfdσ = [ β, 2*σ[2] ]
        end

        # dfdσmax = (βres-βini)/ft*(2*σmax-σ[1]) - βres
        # dfdσmax = (βini-βres)/ft*(σ[1]-σmax)*α*(σmax/ft)^(α-1) - β
        dfdσmax = (βini-βres)/ft*(σ[1]-σmax)*α*(σmax/ft)^(α-1) - β
        m = deriv_σmax_upa(mat, state, up)
        dσmaxdΔλ = m*(norm_r + Δλ*dot(r/norm_r, drdΔλ))
        dfdΔλ = dot(dfdσ, dσdΔλ) + dfdσmax*dσmaxdΔλ
        Δλ = Δλ - f/dfdΔλ
        
        abs(f) < tol && break

        if i == maxits || isnan(Δλ)
            return 0.0, failure("PJoint: Could nof find Δλ.")
        end
    end

    return Δλ, success()
end


function mountD(mat::PJoint, state::PJointState)
    ndim = state.env.ndim
    kn, ks = calc_kn_ks(mat, state)
    α = mat.α
    σmax = calc_σmax(mat, state, state.up)

    De = diagm([kn, ks, ks][1:ndim])

    if state.Δλ == 0.0  # Elastic 
        return De
    elseif σmax == 0.0 && state.w[1] >= 0.0

        Dep = De*1e-4
        return Dep
    else
        # @show "plastic De"
        fc, ft = mat.fc, mat.ft
        βini = 2*ft - fc -2*√(ft^2 - fc*ft)
        βres = mat.γ*βini

        r = potential_derivs(mat, state, state.σ)
        v = yield_derivs(mat, state, state.σ, σmax)
        # dfdσmax = (βres-βini)/ft*(2*σmax-state.σ[1]) - βres
        β = βres + (βini-βres)*(σmax/ft)^mat.α
        dfdσmax = (βini-βres)/ft*(state.σ[1]-σmax)*α*(σmax/ft)^(α-1) - β
        # dfdσmax = -β  # ∂F/∂σmax
        m = deriv_σmax_upa(mat, state, state.up)  # ∂σmax/∂up

        #Dep  = De - De*r*v'*De/(v'*De*r - dfdσmax*m*norm(r))

        if ndim == 3
            den = kn*r[1]*v[1] + ks*r[2]*v[2] + ks*r[3]*v[3] - dfdσmax*m*norm(r)

            Dep = [   kn - kn^2*r[1]*v[1]/den    -kn*ks*r[1]*v[2]/den      -kn*ks*r[1]*v[3]/den
                     -kn*ks*r[2]*v[1]/den         ks - ks^2*r[2]*v[2]/den  -ks^2*r[2]*v[3]/den
                     -kn*ks*r[3]*v[1]/den        -ks^2*r[3]*v[2]/den        ks - ks^2*r[3]*v[3]/den ]
        else
            den = kn*r[1]*v[1] + ks*r[2]*v[2] - dfdσmax*m*norm(r)

            Dep = [   kn - kn^2*r[1]*v[1]/den    -kn*ks*r[1]*v[2]/den      
                     -kn*ks*r[2]*v[1]/den         ks - ks^2*r[2]*v[2]/den  ]
        end

        return Dep
    end
end


function update_state(mat::PJoint, state::PJointState, Δw::Array{Float64,1})

    ndim = state.env.ndim
    σini = copy(state.σ)

    kn, ks = calc_kn_ks(mat, state)
    De = diagm([kn, ks, ks][1:ndim])
    σmax = calc_σmax(mat, state, state.up)  
    # @show σmax

    if isnan(Δw[1]) || isnan(Δw[2])
        alert("PJoint: Invalid value for joint displacement: Δw = $Δw")
    end

    # σ trial and F trial
    σtr  = state.σ + De*Δw

    Ftr  = yield_func(mat, state, σtr, σmax)

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

        # state.σ, state.up = calc_σ_upa(mat, state, σtr)
        if ndim == 3
            if σtr[1] > 0
                state.σ = [σtr[1]/(1 + 2*state.Δλ*kn), σtr[2]/(1 + 2*state.Δλ*ks), σtr[3]/(1 + 2*state.Δλ*ks)]
            else
                state.σ = [σtr[1], σtr[2]/(1 + 2*state.Δλ*ks), σtr[3]/(1 + 2*state.Δλ*ks)]
            end    
        else
            if σtr[1] > 0
                state.σ = [σtr[1]/(1 + 2*state.Δλ*kn), σtr[2]/(1 + 2*state.Δλ*ks)]
            else
                state.σ = [σtr[1], σtr[2]/(1 + 2*state.Δλ*ks)]
            end    
        end
        r = potential_derivs(mat, state, state.σ)
        state.up += state.Δλ*norm(r)

    end
    state.w += Δw
    Δσ = state.σ - σini
    return Δσ, success()
end


function ip_state_vals(mat::PJoint, state::PJointState)
    ndim = state.env.ndim
    if ndim == 3
       return Dict(
          :jw1  => state.w[1],
          :jw2  => state.w[2],
          :jw3  => state.w[3],
          :js1  => state.σ[1],
          :js2  => state.σ[2],
          :js3  => state.σ[3],
          :jup => state.up
          )
    else
        return Dict(
          :jw1  => state.w[1],
          :jw2  => state.w[2],
          :js1  => state.σ[1],
          :js2  => state.σ[2],
          :jup => state.up
          )
    end
end


function output_keys(mat::PJoint)
    return Symbol[:jw1, :js1, :jup]
end