 #This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export TCFJoint

abstract type AbstractTCFJointState<:IpState end

mutable struct TCFJointState<:AbstractTCFJointState
    env::ModelEnv
    σ  ::Array{Float64,1} # stress
    w  ::Array{Float64,1} # relative displacements
    up ::Float64          # effective plastic relative displacement
    Δλ ::Float64          # plastic multiplier
    h  ::Float64          # characteristic length from bulk elements
    w_rate::Float64       # w rate wrt T
    function TCFJointState(env::ModelEnv=ModelEnv())
        this = new(env)
        ndim = env.ndim
        this.σ   = zeros(ndim)
        this.w   = zeros(ndim)
        this.up  = 0.0
        this.Δλ  = 0.0
        this.h   = 0.0
        return this
    end
end

abstract type AbstractTCFJoint<:MatParams end

mutable struct TCFJoint<:AbstractTCFJoint
    E     ::Float64 # Young's modulus
    ν     ::Float64 # Poisson ratio
    ft    ::Float64 # tensile strength
    fc    ::Float64 # compressive strength
    μres  ::Float64 # surface friction after cracking
    ζ     ::Float64 # factor ζ controls the elastic relative displacements (formerly α)
    wc    ::Float64 # critical crack opening
    scurve::String  # softening curve model ("linear" or bilinear" or "hordijk" or "soft")
    α     ::Float64 # curvature coefficient
    τmax0 ::Float64 # initial shear strength limit for σn=0

    function TCFJoint(prms::Dict{Symbol,Float64})
        return  TCFJoint(;prms...)
    end

    function TCFJoint(;E=NaN, nu=NaN, fc=NaN, ft=NaN, mu=NaN, zeta=NaN, alpha=1.0, wc=NaN, GF=NaN, softcurve="hordijk")
        @check GF>0 || wc>0
        @check E>0.0    
        @check 0<=nu<0.5
        @check fc<0.0  
        @check ft>=0.0  
        @check mu>0.0  
        @check zeta>0.0
        @check softcurve in ("hordijk", "soft")
        
        if isnan(wc)
            wc = round(GF/(0.1947019536*ft), sigdigits=5)  
        end
        α = alpha
        a = (2*α*fc - 2*ft - fc + √(fc^2*(2*α-1)^2 - 4*fc*ft + 4*ft^2)) / (4*α-4)
        b = √((2*a-fc)*(ft-a)/(2*α))
        τmax0 = b/((ft-a)/ft)^α

        this = new(E, nu, ft, fc, mu, zeta, wc, softcurve, α, τmax0)
        return this
    end
end


function paramsdict(matparams::AbstractTCFJoint)
    matparams = OrderedDict( string(field)=> getfield(matparams, field) for field in fieldnames(typeof(matparams)) )

    matparams.scurve in ("hordijk", "soft") && ( matparams["GF"] = 0.1943*matparams.ft*matparams.wc )
    return matparams
end


# Returns the element type that works with this material model
matching_elem_type(::AbstractTCFJoint) = MechJointElem

# Type of corresponding state structure
ip_state_type(::MechJointElem, ::AbstractTCFJoint) = TCFJointState


function yield_func(matparams::AbstractTCFJoint, state::AbstractTCFJointState, σ::Array{Float64,1}, up::Float64)
    σmax = calc_σmax(matparams, state, up)
    τmax = calc_τmax(matparams, state, up)
    μ    = calc_μ(matparams, state, up)
    tmp  = (σmax-σ[1])/matparams.ft

    if state.env.ndim == 3
        τnorm = sqrt(σ[2]^2 + σ[3]^2)
    else
        τnorm = abs(σ[2])
    end

    if σmax>σ[1]
        return τnorm - τmax*tmp^matparams.α  - (σmax - σ[1])*μ
    else
        return τnorm - (σmax - σ[1])*μ
    end
end


function yield_derivs(matparams::AbstractTCFJoint, state::AbstractTCFJointState, σ::Array{Float64,1}, up::Float64)
    ft    = matparams.ft
    α     = matparams.α
    σmax  = calc_σmax(matparams, state, up)
    τmax  = calc_τmax(matparams, state, up)
    μ     = calc_μ(matparams, state, up)
    tmp   = (σmax-σ[1])/ft

    if σmax>σ[1]
        dfdσn   = μ + α*τmax/ft*tmp^(α-1)
        dfdσmax = -μ - α*τmax/ft*tmp^(α-1)
        dfdτmax = -tmp^α
    else
        dfdσn   = μ
        dfdσmax = -μ
        dfdτmax = 0.0
    end
    
    dfdμ = -(σmax-σ[1])

    if state.env.ndim == 3
        τnorm = sqrt(σ[2]^2 + σ[3]^2)
        dfdσ = [ dfdσn, σ[2]/τnorm, σ[3]/τnorm]
    else
        dfdσ = [ dfdσn, sign(σ[2]) ]
    end

    return dfdσ, dfdσmax, dfdτmax, dfdμ
end


function potential_derivs(matparams::AbstractTCFJoint, state::AbstractTCFJointState, σ::Array{Float64,1})
    ndim = state.env.ndim
    if ndim == 3
        if σ[1] > 0.0 
            # G1:
            r = Float64[ 2.0*σ[1], 2.0*σ[2], 2.0*σ[3]]
        else
            # G2:
            r = Float64[ 0.0, 2.0*σ[2], 2.0*σ[3] ]
        end
        if r[1]==r[2]==r[3]==0.0
            r = Float64[ 1.0, 0.0, 0.0 ] # important
        end
    else
        if σ[1] > 0.0 
            # G1:
            r = Float64[ 2*σ[1], 2*σ[2]]
        else
            # G2:
            r = Float64[ 0.0, 2*σ[2] ]
        end
        if r[1]==r[2]==0.0
            r = Float64[ 1.0, 0.0 ] # important
        end
    end
    return r
end


function calc_σmax(matparams::AbstractTCFJoint, state::AbstractTCFJointState, up::Float64)
    if matparams.scurve == "hordijk"
        if up < matparams.wc
            z = (1 + 27*(up/matparams.wc)^3)*exp(-6.93*up/matparams.wc) - 28*(up/matparams.wc)*exp(-6.93)
        else
            z = 0.0
        end
    elseif matparams.scurve == "soft"
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
    end

    return matparams.ft*z
end


function deriv_σmax_up(matparams::AbstractTCFJoint, state::AbstractTCFJointState, up::Float64)
    x = up/matparams.wc
    if matparams.scurve == "hordijk"
        if x<1
            dzdx = 81*x^2*exp(-6.93*x) - 6.93*(1 + 27*x^3)*exp(-6.93*x) - 0.02738402432
        else
            dzdx = 0.0
        end
    elseif matparams.scurve == "soft"
        b = 1.30837
        α = 0.55

        if up == 0.0
            dzdx = 0.0
        elseif up < matparams.wc
            dzdx = -α*log(b)*b^(1-x^-α)/x^(α+1)
        else
            dzdx = 0.0
        end
    end

    return matparams.ft*dzdx/matparams.wc
end


function calc_τmax(matparams::AbstractTCFJoint, state::AbstractTCFJointState, up::Float64)
    xp = 0.13
    y0 = 0.8
    y0 = 0.7

    y0 = 0.7
    xp = 0.08
    
    y0 = 0.9
    xp = 0.01


    b  = 1.3
    b  = exp(1)
    α2 = 0.5

    α2 = 0.55
    b = 1.30837
    
    α2 = 0.5
    b = 2.0


    x = up/matparams.wc

    if x<=xp
        α1 = 0.4
        z  = y0 + (1-y0)*(x/xp)^α1
    elseif x<1
        x̄ = (x-xp)/(1-xp)
        z = 1 - b^(1-x̄^-α2)
    else
        z = 0.0
    end

    return matparams.τmax0*z
    
end


function deriv_τmax_up(matparams::AbstractTCFJoint, state::AbstractTCFJointState, up::Float64)
    xp = 0.13
    y0 = 0.8
    y0 = 0.7

    y0 = 0.7
    xp = 0.08

    y0 = 0.9
    xp = 0.01

    # b  = 1.3
    b  = exp(1)
    α2 = 0.5

    α2 = 0.55
    b = 1.30837

    α2 = 0.5
    b = 2.0

    x = up/matparams.wc
    
    if x<=xp
        xmin = 1e-4
        x    = max(x, xmin)
        α1   = 0.4
        dzdx = α1*(1-y0)*(x/xp)^(α1-1)/xp
    elseif x<1
        x̄ = (x-xp)/(1-xp)
        dzdx̄ = -α2*log(b)*b^(1-x̄^-α2)/x̄^(α2+1)
        dx̄dx = 1/(1-xp)
        dzdx = dzdx̄*dx̄dx
    else
        dzdx = 0.0
    end

    return matparams.τmax0*dzdx/matparams.wc
end


function calc_μ(matparams::AbstractTCFJoint, state::AbstractTCFJointState, up::Float64)
    x = up/matparams.wc
    if x==0
        z = 1e-10
    elseif x<1
        z = x
    else
        z = 1.0
    end
    return matparams.μres*z
end


function deriv_μ_up(matparams::AbstractTCFJoint, state::AbstractTCFJointState, up::Float64)
    x = up/matparams.wc
    if x<1.0
        dzdx = 1.0
    else
        dzdx = 0.0
    end
    return matparams.μres*dzdx/matparams.wc
end


function calc_kn_ks(matparams::AbstractTCFJoint, state::AbstractTCFJointState)
    kn = matparams.E*matparams.ζ/state.h
    G  = matparams.E/(2.0*(1.0+matparams.ν))
    ks = G*matparams.ζ/state.h

    return kn, ks
end


function mountD(matparams::AbstractTCFJoint, state::AbstractTCFJointState)

    ndim   = state.env.ndim
    kn, ks = calc_kn_ks(matparams, state)
    σmax   = calc_σmax(matparams, state, state.up)

    De = diagm([kn, ks, ks][1:ndim])

    if state.Δλ == 0.0  # Elastic 
        return De
    elseif σmax == 0.0 && state.w[1] >= 0.0
        Dep = De*1e-4
        # Dep = De*1e-3
        return Dep
    else
        up = state.up
        v, dfdσmax, dfdτmax, dfdμ = yield_derivs(matparams, state, state.σ, up)
        r = potential_derivs(matparams, state, state.σ)
        dσmaxdup = deriv_σmax_up(matparams, state, up)
        dτmaxdup = deriv_τmax_up(matparams, state, up)
        dμdup    = deriv_μ_up(matparams, state, up)
        dfdup    = dfdσmax*dσmaxdup + dfdτmax*dτmaxdup + dfdμ*dμdup

        if ndim == 3
            den = kn*r[1]*v[1] + ks*r[2]*v[2] + ks*r[3]*v[3] - dfdup*norm(r)

            Dep = [   kn - kn^2*r[1]*v[1]/den    -kn*ks*r[1]*v[2]/den      -kn*ks*r[1]*v[3]/den
                     -kn*ks*r[2]*v[1]/den         ks - ks^2*r[2]*v[2]/den  -ks^2*r[2]*v[3]/den
                     -kn*ks*r[3]*v[1]/den        -ks^2*r[3]*v[2]/den        ks - ks^2*r[3]*v[3]/den ]
        else
            den = kn*r[1]*v[1] + ks*r[2]*v[2] - dfdup*norm(r)

            Dep = [   kn - kn^2*r[1]*v[1]/den    -kn*ks*r[1]*v[2]/den      
                     -kn*ks*r[2]*v[1]/den         ks - ks^2*r[2]*v[2]/den  ]
        end

        return Dep
    end
end


function calc_σ_up_Δλ(matparams::AbstractTCFJoint, state::AbstractTCFJointState, σtr::Array{Float64,1})

    ndim = state.env.ndim
    Δλ   = 0.0
    up   = 0.0
    σ    = zeros(ndim)
    σ0   = zeros(ndim)

    tol    = 1e-4 # better
    tol    = 1e-5 
    tol    = 1e-6

    tol = 10^-(10-log10(matparams.ft))

    maxits = 50
    for i in 1:maxits
        kn, ks = calc_kn_ks(matparams, state)

        # quantities at n+1
        if ndim == 3
            if σtr[1]>0
                 σ     = Float64[ σtr[1]/(1+2*Δλ*kn),  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
                 dσdΔλ = Float64[ -2*kn*σtr[1]/(1+2*Δλ*kn)^2,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2,  -2*ks*σtr[3]/(1+2*Δλ*ks)^2 ]
            else
                 σ     = Float64[ σtr[1],  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
                 dσdΔλ = Float64[ 0,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2,  -2*ks*σtr[3]/(1+2*Δλ*ks)^2 ]
            end
        else
            if σtr[1]>0
                 σ     = Float64[ σtr[1]/(1+2*Δλ*kn),  σtr[2]/(1+2*Δλ*ks) ]
                 dσdΔλ = Float64[ -2*kn*σtr[1]/(1+2*Δλ*kn)^2,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2 ]
            else
                 σ     = Float64[ σtr[1],  σtr[2]/(1+2*Δλ*ks) ]
                 dσdΔλ = Float64[ 0,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2 ]
             end
        end

        drdΔλ = 2*dσdΔλ

        r  = potential_derivs(matparams, state, σ)
        nr = norm(r)
        up = state.up + Δλ*nr
       
        f  = yield_func(matparams, state, σ, up)
        dfdσ, dfdσmax, dfdτmax, dfdμ = yield_derivs(matparams, state, σ, up)
        
        dσmaxdup  = deriv_σmax_up(matparams, state, up)
        dτmaxdup  = deriv_τmax_up(matparams, state, up)
        dμdup     = deriv_μ_up(matparams, state, up)
        dupdΔλ    = nr + Δλ*dot(r/nr, drdΔλ)

        dfdup = dfdσmax*dσmaxdup + dfdτmax*dτmaxdup + dfdμ*dμdup
        dfdΔλ = dot(dfdσ, dσdΔλ) + dfdup*dupdΔλ
        Δλ = Δλ - f/dfdΔλ

        if Δλ<=0 || isnan(Δλ) || i==maxits
            # switch to bissection method
            return calc_σ_up_Δλ_bissection(matparams, state, σtr)
        end

        if maximum(abs, σ-σ0) <= tol
        # if abs(f) <= tol
            break
        end
        σ0 .= σ
    end

    return σ, up, Δλ, success()   
end


function yield_func_from_Δλ(matparams::AbstractTCFJoint, state::AbstractTCFJointState, σtr::Array{Float64,1}, Δλ::Float64)
    ndim    = state.env.ndim
    kn, ks  = calc_kn_ks(matparams, state)

    # quantities at n+1
    if ndim == 3
        if σtr[1]>0
            σ = Float64[ σtr[1]/(1+2*Δλ*kn),  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
        else
            σ = Float64[ σtr[1],  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
        end
    else
        if σtr[1]>0
            σ = Float64[ σtr[1]/(1+2*Δλ*kn),  σtr[2]/(1+2*Δλ*ks) ]
        else
            σ = Float64[ σtr[1],  σtr[2]/(1+2*Δλ*ks) ]
        end
    end

    r  = potential_derivs(matparams, state, σ)
    nr = norm(r)
    up = state.up + Δλ*nr
    f  = yield_func(matparams, state, σ, up)
    return f
end


function calc_σ_up_Δλ_bissection(matparams::AbstractTCFJoint, state::AbstractTCFJointState, σtr::Array{Float64,1})
    ndim    = state.env.ndim
    kn, ks  = calc_kn_ks(matparams, state)
    De      = diagm([kn, ks, ks][1:ndim])
    r       = potential_derivs(matparams, state, state.σ)

    # Δλ estimative
    Δλ0 = norm(σtr-state.σ)/norm(De*r)
    
    # find initial interval
    a  = 0.0
    b  = Δλ0
    fa = yield_func_from_Δλ(matparams, state, σtr, a)
    fb = yield_func_from_Δλ(matparams, state, σtr, b)

    # search for a valid interval
    if fa*fb>0
        δ = Δλ0
        maxits = 100
        for i in 1:maxits
            b += δ
            fb = yield_func_from_Δλ(matparams, state, σtr, b)
            if fa*fb<0.0
                break
            end
            δ *= 1.5

            if i==maxits
                return 0.0, state.σ, 0.0, failure("TCFJoint: could not find Δλ (could not find interval for bissection)")
            end
        end
    end

    # bissection method
    Δλ = 0.0
    f  = 0.0
    up = 0.0
    σ0 = zeros(ndim) # previous value
    σ  = zeros(ndim)

    tol = 10^-(10-log10(matparams.ft))

    # tol = 1e-4
    # tol = 1e-5
    # tol = 1e-6
    # tol = 1e-8

    maxits = 100
    for i in 1:maxits
        Δλ = (a+b)/2

        f = yield_func_from_Δλ(matparams, state, σtr, Δλ)

        # quantities at n+1
        if ndim == 3
            if σtr[1]>0
                σ = Float64[ σtr[1]/(1+2*Δλ*kn),  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
            else
                σ = Float64[ σtr[1],  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
            end
        else
            if σtr[1]>0
                σ = Float64[ σtr[1]/(1+2*Δλ*kn),  σtr[2]/(1+2*Δλ*ks) ]
            else
                σ = Float64[ σtr[1],  σtr[2]/(1+2*Δλ*ks) ]
            end
        end

        r  = potential_derivs(matparams, state, σ)
        nr = norm(r)
        up = state.up + Δλ*nr
        f  = yield_func(matparams, state, σ, up)

        # @show f
        # @show Δλ
        # error()

        if fa*f<0
            b = Δλ
        else
            a  = Δλ
            fa = f
        end

        if maximum(abs, σ-σ0) <= tol
            break
        end
        σ0 .= σ

        i==maxits && return 0.0, state.σ, 0.0, failure("TCFJoint: could not find Δλ with NR/bissection (maxits reached, f=$f)")
    end

    return σ, up, Δλ, success()   
end


function find_σ_up_Δλ(matparams::AbstractTCFJoint, state::AbstractTCFJointState, σtr::Array{Float64,1})
    ndim    = state.env.ndim
    kn, ks  = calc_kn_ks(matparams, state)
    De      = diagm([kn, ks, ks][1:ndim])
    r       = potential_derivs(matparams, state, state.σ)

    # Δλ estimative
    Δλ0 = norm(σtr-state.σ)/norm(De*r)
    
    # find initial interval
    a  = 0.0
    b  = Δλ0
    fa = yield_func_from_Δλ(matparams, state, σtr, a)
    fb = yield_func_from_Δλ(matparams, state, σtr, b)

    # search for a valid interval
    if fa*fb>0
        δ = Δλ0
        maxits = 100
        for i in 1:maxits
            b += δ
            fb = yield_func_from_Δλ(matparams, state, σtr, b)
            fa*fb<0.0 && break
            δ *= 1.5
            i==maxits && return 0.0, state.σ, 0.0, failure("TCFJoint: could not find Δλ (could not find interval for bissection)")
        end
    end

    tol   = Δλ0*10^-3
    f(Δλ) = yield_func_from_Δλ(matparams, state, σtr, Δλ)
    Δλ, status = findroot(f, a, b, tol)
    failed(status) && return state.σ, 0.0, 0.0, status

    # quantities at n+1
    if ndim == 3
        if σtr[1]>0
            σ = Float64[ σtr[1]/(1+2*Δλ*kn),  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
        else
            σ = Float64[ σtr[1],  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
        end
    else
        if σtr[1]>0
            σ = Float64[ σtr[1]/(1+2*Δλ*kn),  σtr[2]/(1+2*Δλ*ks) ]
        else
            σ = Float64[ σtr[1],  σtr[2]/(1+2*Δλ*ks) ]
        end
    end

    r  = potential_derivs(matparams, state, σ)
    nr = norm(r)
    up = state.up + Δλ*nr

    return σ, up, Δλ, success()   
end


function update_state(matparams::AbstractTCFJoint, state::AbstractTCFJointState, Δw::Array{Float64,1})
    ndim = state.env.ndim
    σini = copy(state.σ)

    kn, ks = calc_kn_ks(matparams, state)
    De = diagm([kn, ks, ks][1:ndim])

    if isnan(Δw[1]) || isnan(Δw[2])
        alert("AbstractTCFJoint: Invalid value for joint displacement: Δw = $Δw")
    end

    # σ trial and F trial
    σtr  = state.σ + De*Δw

    Ftr  = yield_func(matparams, state, σtr, state.up)

    # Elastic and EP integration
    if Ftr <= 0.0
        state.Δλ  = 0.0
        state.σ  .= σtr
    elseif state.up>=matparams.wc && σtr[1]>0
        if ndim==3
            Δup = norm([ σtr[1]/kn, σtr[2]/ks, σtr[3]/ks ])
        else
            Δup = norm([ σtr[1]/kn, σtr[2]/ks ])
        end
        state.up += Δup
        state.σ  .= 0.0
        state.Δλ  = 1.0
    else
        # Plastic increment
        # σ, up, Δλ, status = calc_σ_up_Δλ(matparams, state, σtr)
        σ, up, Δλ, status = find_σ_up_Δλ(matparams, state, σtr)
        failed(status) && return state.σ, status

        state.σ, state.up, state.Δλ = σ, up, Δλ
    end
    state.w += Δw
    Δσ = state.σ - σini
    return Δσ, success()
end


function ip_state_vals(matparams::AbstractTCFJoint, state::AbstractTCFJointState)
    ndim = state.env.ndim
    if ndim == 3
       return Dict(
          :jw1 => state.w[1],
          :jw2 => state.w[2],
          :jw3 => state.w[3],
          :js1 => state.σ[1],
          :js2 => state.σ[2],
          :js3 => state.σ[3],
          :jup => state.up
          )
    else
        return Dict(
          :jw1 => state.w[1],
          :jw2 => state.w[2],
          :js1 => state.σ[1],
          :js2 => state.σ[2],
          :jup => state.up
          )
    end
end


function output_keys(matparams::AbstractTCFJoint)
    return Symbol[:jw1, :js1, :jup]
end