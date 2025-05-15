 #This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export MMCJoint

mutable struct MMCJointState<:IpState
    ctx::Context
    σ  ::Array{Float64,1} # stress
    w  ::Array{Float64,1} # relative displacements
    up::Float64          # effective plastic relative displacement
    Δλ ::Float64          # plastic multiplier
    h  ::Float64          # characteristic length from bulk elements
    function MMCJointState(ctx::Context)
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

mutable struct MMCJoint<:Material
    E ::Float64      # Young's modulus
    ν ::Float64      # Poisson ratio
    ζ ::Float64      # factor ζ controls the elastic relative displacements (formerly α)
    ft::Float64      # tensile strength (internal variable)
    fc::Float64      # compressive strength (internal variable)
    α ::Float64      # exponent related to the yield function opening
    β ::Float64      # coefficient used to match the compression strenght
    wc::Float64      # critical crack opening
    ws::Float64      # openning at inflection (where the curve slope changes)
    softcurve::String # softening curve model ("linear" or bilinear" or "hordijk")

    function MMCJoint(prms::Dict{Symbol,Float64})
        return  MMCJoint(;prms...)
    end

    function MMCJoint(;E=NaN, nu=NaN, fc=NaN, ft=NaN, zeta=NaN, alpha=1.0, wc=NaN, ws=NaN, GF=NaN, Gf=NaN, softcurve="bilinear")

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
        @check 0.5<=alpha
        zeta>0     || error("Invalid value for zeta: $zeta")
        wc>0       || error("Invalid value for wc: $wc")
        (isnan(ws) || ws>0) || error("Invalid value for ws: $ws")
        softcurve in ("linear", "bilinear", "hordijk") || error("Invalid softcurve: softcurve must to be linear or bilinear or hordijk")

        # (a,b) is the point for simple compression failure
        a = (2*alpha*ft + alpha*fc-fc-√( (alpha*fc-fc+2*alpha*ft)^2 + 4*alpha*fc*ft*(1-2*alpha) )) / (4*alpha-2)
        b = √(alpha*(2*a-fc)*(ft-a))
        beta = (ft-a)/(b^2/ft^2)^alpha

        this = new(E, nu, zeta, ft, fc, alpha, beta, wc, ws, softcurve)
        return this
    end
end


# Type of corresponding state structure
compat_state_type(::Type{MMCJoint}) = MMCJointState

# Element types that work with this material
compat_elem_types(::Type{MMCJoint}) = (MechJoint,)


function yield_func(mat::MMCJoint, state::MMCJointState, σ::Array{Float64,1}, σmax::Float64)
    ft, α, β = mat.ft, mat.α, mat.β

    if state.ctx.ndim == 3
        return σ[1] - σmax + β*((σ[2]^2 + σ[3]^2)/ft^2)^α
    else
        return σ[1] - σmax + β*(σ[2]^2/ft^2)^α
    end
end


function yield_derivs(mat::MMCJoint, state::MMCJointState, σ::Array{Float64,1})
    ft, α, β = mat.ft, mat.α, mat.β

    if state.ctx.ndim == 3
        tmp = 2*α*β/ft^2*((σ[2]^2+σ[3]^2)/ft^2)^(α-1)
        return [ 1 , σ[2]*tmp, σ[3]*tmp ]
    else
        tmp = 2*α*β/ft^2*(σ[2]^2/ft^2)^(α-1)
        return [ 1 , σ[2]*tmp ]
    end
end


function potential_derivs(mat::MMCJoint, state::MMCJointState, σ::Array{Float64,1})
    ndim = state.ctx.ndim
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


function calc_σmax(mat::MMCJoint, state::MMCJointState, up::Float64)
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


function deriv_σmax_upa(mat::MMCJoint, state::MMCJointState, up::Float64)
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


function calc_kn_ks(mat::MMCJoint, state::MMCJointState)
    kn = mat.E*mat.ζ/state.h
    G  = mat.E/(2.0*(1.0+mat.ν))
    ks = G*mat.ζ/state.h

    return kn, ks
end


function calc_Δλ(mat::MMCJoint, state::MMCJointState, σtr::Array{Float64,1})
    ndim = state.ctx.ndim
    maxits = 20
    Δλ     = 0.0
    f      = 0.0
    up    = 0.0
    tol    = 1e-4
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

        f = yield_func(mat, state, σ, σmax)
        dfdσ = yield_derivs(mat, state, σ)

        m = deriv_σmax_upa(mat, state, up)
        dσmaxdΔλ = m*(norm_r + Δλ*dot(r/norm_r, drdΔλ))
        dfdσmax = -1
        dfdΔλ = dot(dfdσ, dσdΔλ) + dfdσmax*dσmaxdΔλ
        Δλ = Δλ - f/dfdΔλ
        # @show f, Δλ
        abs(f) < tol && break

        if i == maxits || isnan(Δλ)
            # warn("""MMCJoint: Could not find Δλ. This may happen when the system
            # becomes hypostatic and thus the global stiffness matrix is nearly singular.
            # Increasing the mesh refinement may result in a nonsingular matrix.
            # """)
            # warn("iterations=$i Δλ=$Δλ")
            return 0.0, failure("MMCJoint: Could nof find Δλ.")
        end
    end
    # @show nits
    # @show f, Δλ

    return Δλ, success()
end


function calc_σ_upa(mat::MMCJoint, state::MMCJointState, σtr::Array{Float64,1})
    ndim = state.ctx.ndim
    kn, ks = calc_kn_ks(mat, state)

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
    return state.σ, state.up
end


function calcD(mat::MMCJoint, state::MMCJointState)
    ndim = state.ctx.ndim
    kn, ks = calc_kn_ks(mat, state)
    σmax = calc_σmax(mat, state, state.up)

    De = diagm([kn, ks, ks][1:ndim])
    # @show σmax
    # @show state.up
    # @show state.up > mat.wc

    if state.Δλ == 0.0  # Elastic 
        # @show "ELASTIC De"
        return De
    # elseif σmax == 0.0 # ?????????????? and w1>0 ???
    elseif σmax == 0.0 && state.w[1]>0
        # @show "smax=0 De"

        Dep = De*1e-4
        # Dep  = diagm(ones(ndim))
        # Dep[1,1] *= -1
        return Dep
    else
        # @show "plastic De"
        r = potential_derivs(mat, state, state.σ)
        v = yield_derivs(mat, state, state.σ)
        y = -1  # ∂F/∂σmax
        m = deriv_σmax_upa(mat, state, state.up)  # ∂σmax/∂up

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


        # if σmax == 0.0 
        #     # @show Dep*den
        #     # Dep += 1e-8*De
        #     # @show De
        #     @show kn
        #     @show ks
        #     @show r
        #     @show v
        #     @show y
        #     @show m
        #     @show den
        #     @show Dep
        #     error()
        # end

        return Dep
    end
end


function update_state!(mat::MMCJoint, state::MMCJointState, Δw::Array{Float64,1})

    ndim = state.ctx.ndim
    σini = copy(state.σ)

    kn, ks = calc_kn_ks(mat, state)
    De = diagm([kn, ks, ks][1:ndim])
    σmax = calc_σmax(mat, state, state.up)  
    # @show σmax

    if isnan(Δw[1]) || isnan(Δw[2])
        alert("MMCJoint: Invalid value for joint displacement: Δw = $Δw")
    end

    # σ trial and F trial
    σtr  = state.σ + De*Δw

    Ftr  = yield_func(mat, state, σtr, σmax)
    # @show "stress update"
    # @show σmax
    # @show state.up
    # @show Δw
    # @show state.σ
    # @show σtr
    # @show Ftr
    # @show state.w[1] 

    # Elastic and EP integration
    if σmax == 0.0 && state.w[1] >= 0.0
        # @show "smax=0 up"
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
        # @show "ELASTIC up"
        # Pure elastic increment
        state.Δλ = 0.0
        state.σ  = copy(σtr) 

    else
        # @show "plastic up"

        # Plastic increment
        state.Δλ, status = calc_Δλ(mat, state, σtr)
        failed(status) && return state.σ, status

        state.σ, state.up = calc_σ_upa(mat, state, σtr)

        # @show state.Δλ
                      
        # Return to surface:
        # F  = yield_func(mat, state, state.σ)   
        # F > 1e-2 && alert("MMCJoint: Yield function value ($F) outside tolerance")

    end
    state.w += Δw
    Δσ = state.σ - σini
    return Δσ, success()
end


function ip_state_vals(mat::MMCJoint, state::MMCJointState)
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


function output_keys(mat::MMCJoint)
    return Symbol[:jw, :jσn, :jup]
end