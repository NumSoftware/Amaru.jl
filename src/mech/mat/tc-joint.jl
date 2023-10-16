 #This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export TCJoint

mutable struct TCJointState<:IpState
    env::ModelEnv
    σ  ::Array{Float64,1} # stress
    w  ::Array{Float64,1} # relative displacements
    up ::Float64          # effective plastic relative displacement
    Δλ ::Float64          # plastic multiplier
    h  ::Float64          # characteristic length from bulk elements
    function TCJointState(env::ModelEnv)
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

mutable struct TCJoint<:Material
    E ::Float64       # Young's modulus
    ν ::Float64       # Poisson ratio
    ft::Float64       # tensile strength
    fc::Float64       # compressive strength
    ζ ::Float64       # factor ζ controls the elastic relative displacements (formerly α)
    wc::Float64       # critical crack opening
    ws::Float64       # openning at inflection (where the curve slope changes)
    softcurve::String # softening curve model ("linear" or bilinear" or "hordijk" or "soft")
    α::Float64        # curvature coefficient
    γ::Float64        # factor for βres
    θ::Float64        # fator for the surface reduction speed
    βini::Float64     # initial curvature size

    function TCJoint(prms::Dict{Symbol,Float64})
        return  TCJoint(;prms...)
    end

    function TCJoint(;E=NaN, nu=NaN, fc=NaN, ft=NaN, zeta=NaN, alpha=1.5, gamma=0.1, theta=1.5, wc=NaN, ws=NaN, GF=NaN, Gf=NaN, softcurve="hordijk")
        @check GF>0 || Gf>0 || wc>0
        @check E>0.0    
        @check 0<=nu<0.5
        @check fc<0  
        @check ft>=0  
        @check zeta>0
        @check softcurve in ("linear", "bilinear", "hordijk", "soft")
        
        if isnan(wc)
            if softcurve == "linear"
                wc = round(2*GF/ft, sigdigits=5)
            elseif softcurve == "bilinear"
                if isnan(Gf)
                    wc = round(5*GF/ft, sigdigits=5)
                    ws = round(wc*0.15, sigdigits=5)
                else
                    wc = round((8*GF- 6*Gf)/ft, sigdigits=5)
                    ws = round(1.5*Gf/ft, sigdigits=5)
                end
            elseif softcurve=="hordijk"
                wc = round(GF/(0.1947019536*ft), sigdigits=5)  
            elseif softcurve=="soft"
                wc = round(GF/(0.1947019536*ft), sigdigits=5)  
                # wc = round(6*GF/ft, sigdigits=5)  
                # wc = round(5.14*GF/ft, sigdigits=5)  
            end    
        end
        
        @check isnan(ws) || ws>0

        a = (2*alpha*ft + alpha*fc - fc - √(alpha^2*fc^2 - 4*alpha^2*fc*ft + 4*alpha^2*ft^2 - 2*alpha*fc^2 + fc^2)) / (4*alpha-2)
        b = √(alpha*(2*a-fc)*(ft-a))
        βini = (b^2/ft^2)^alpha/(ft-a)

        this = new(E, nu, ft, fc, zeta, wc, ws, softcurve, alpha, gamma, theta, βini)
        return this
    end
end


function paramsdict(mat::TCJoint)
    mat = OrderedDict( string(field)=> getfield(mat, field) for field in fieldnames(typeof(mat)) )

    mat.softcurve in ("hordijk", "soft") && ( mat["GF"] = 0.1943*mat.ft*mat.wc )
    return mat
end


# Type of corresponding state structure
compat_state_type(::Type{TCJoint}, ::Type{MechJoint}, env::ModelEnv) = TCJointState

# Element types that work with this material
# compat_elem_types(::Type{TCJoint}) = (MechJoint,)


function beta(mat::TCJoint, σmax::Float64)
    βini = mat.βini
    βres = mat.γ*βini
    return βres + (βini-βres)*(σmax/mat.ft)^mat.α
end


function yield_func(mat::TCJoint, state::TCJointState, σ::Array{Float64,1}, σmax::Float64)
    α  = mat.α
    β = beta(mat, σmax)
    ft = mat.ft
    if state.env.ndim == 3
        return β*(σ[1] - σmax) + ((σ[2]^2 + σ[3]^2)/ft^2)^α
    else
        return β*(σ[1] - σmax) + (σ[2]^2/ft^2)^α
    end
end


function yield_derivs(mat::TCJoint, state::TCJointState, σ::Array{Float64,1}, σmax::Float64)
    α = mat.α
    β = beta(mat, σmax)
    ft = mat.ft

    if state.env.ndim == 3
        tmp = 2*α/ft^2*((σ[2]^2+σ[3]^2)/ft^2)^(α-1)
        σ[2]==σ[3]==0.0 && (tmp=0)
        return [ β , σ[2]*tmp, σ[3]*tmp ]
    else
        tmp = 2*α/ft^2*(σ[2]^2/ft^2)^(α-1)
        σ[2]==0.0 && (tmp=0)
        return [ β , σ[2]*tmp ]
    end
end


function potential_derivs(mat::TCJoint, state::TCJointState, σ::Array{Float64,1})
    ndim = state.env.ndim
    if ndim == 3
        if σ[1] > 0.0 
            # G1:
            r = [ 2.0*σ[1], 2.0*σ[2], 2.0*σ[3]]
        else
            # G2:
            r = [ 0.0, 2.0*σ[2], 2.0*σ[3] ]
        end
        if r[1]==r[2]==r[3]==0.0
            r = [ 1.0, 0.0, 0.0]
        end
    else
        if σ[1] > 0.0 
            # G1:
            r = [ 2*σ[1], 2*σ[2]]
        else
            # G2:
            r = [ 0.0, 2*σ[2] ]
        end
        if r[1]==r[2]==0.0
            r = [ 1.0, 0.0]
        end
    end
    return r
end


function calc_σmax(mat::TCJoint, state::TCJointState, up::Float64)
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
            z = (1 + 27*(up/mat.wc)^3)*exp(-6.93*up/mat.wc) - 28*(up/mat.wc)*exp(-6.93)
        else
            z = 0.0
        end
        σmax = z*mat.ft
    elseif mat.softcurve == "soft"
        # c = 1.5
        # a = 1.19311
        # c = 0.0
        # a = 1.09913
        # c = 2.0
        # a = 1.24787
        # a = 1.10366
        # a = 1.07946
        m = 0.55
        a = 1.30837
        if up == 0.0
            z = 1.0
        elseif 0.0 < up < mat.wc
            x = up/mat.wc
            # z = (1.0 - x)^c*(1.0 - a^(1.0 - 1.0/x))
            z = 1.0 - a^(1.0 - 1.0/x^m)
        else
            z = 0.0
        end
        σmax = z*mat.ft
    end

    # @show mat.softcurve

    return σmax
end


function deriv_σmax_upa(mat::TCJoint, state::TCJointState, up::Float64)
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
            dz = ((81*up^2*exp(-6.93*up/mat.wc)/mat.wc^3) - (6.93*(1 + 27*up^3/mat.wc^3)*exp(-6.93*up/mat.wc)/mat.wc) - 0.02738402432/mat.wc)
        else
            dz = 0.0
        end
        dσmax = dz*mat.ft 
    elseif mat.softcurve == "soft"
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

    # @show mat.softcurve

    return dσmax
end


function calc_kn_ks(mat::TCJoint, state::TCJointState)
    kn = mat.E*mat.ζ/state.h
    G  = mat.E/(2*(1+mat.ν))
    ks = G*mat.ζ/state.h

    return kn, ks
end


function consistentD(mat::TCJoint, state::TCJointState)
    # numerical approximation
    # seems not to work under compressive loads

    ndim = state.env.ndim
    σmax = calc_σmax(mat, state, state.up)

    if state.Δλ == 0.0
        kn, ks = calc_kn_ks(mat, state)
        De = diagm([kn, ks, ks][1:ndim])
        return De
    # elseif σmax == 0.0 && state.w[1] >= 0.0
    #     kn, ks = calc_kn_ks(mat, state)
    #     De = diagm([kn, ks, ks][1:ndim])
    #     Dep = De*1e-4
    #     # Dep = De*1e-3
    #     return Dep
    end

    Dep = zeros(ndim, ndim)
    V = zeros(ndim)
    h = √eps()
    h = eps()^(1/3)

    # iteration for all w components
    for j in 1:ndim
        statej = copy(state)
        V[j] = 1.0
        Δw = h*V
        Δσ, succeeded = update_state!(mat, statej, Δw)
        Dep[:,j] .= Δσ./h
        V[j] = 0.0
    end

    return Dep
end


function calcD(mat::TCJoint, state::TCJointState)
    # return consistentD(mat, state)

    ndim = state.env.ndim
    kn, ks = calc_kn_ks(mat, state)
    θ = mat.θ
    σmax = calc_σmax(mat, state, state.up)

    De = diagm([kn, ks, ks][1:ndim])

    if state.Δλ == 0.0  # Elastic 
        # @show "Elastic"
        return De
    elseif σmax == 0.0 && state.w[1] >= 0.0
        # @show "Plast"
        Dep = De*1e-4
        # Dep = De*1e-3
        return Dep
    else
        # @show "Elastic Pla"

        ft = mat.ft
        βini = mat.βini
        βres = mat.γ*βini
        β = beta(mat, σmax)
        dfdσmax = (βini-βres)/ft*(state.σ[1]-σmax)*θ*(σmax/ft)^(θ-1) - β

        r = potential_derivs(mat, state, state.σ)
        v = yield_derivs(mat, state, state.σ, σmax)
        m = deriv_σmax_upa(mat, state, state.up)  # ∂σmax/∂up

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


function calc_σ_up_Δλ(mat::TCJoint, state::TCJointState, σtr::Array{Float64,1})
    ndim = state.env.ndim
    Δλ   = 0.0
    up   = 0.0
    σ    = zeros(ndim)
    σ0   = zeros(ndim)
    βini = mat.βini
    βres = mat.γ*βini
    θ    = mat.θ
    ft   = mat.ft
    
    tol    = 1e-6
    maxits = 50
    for i in 1:maxits
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
        β      = beta(mat, σmax)

        f    = yield_func(mat, state, σ, σmax)
        dfdσ = yield_derivs(mat, state, σ, σmax)

        dfdσmax = (βini-βres)/ft*(σ[1]-σmax)*θ*(σmax/ft)^(θ-1) - β
        m = deriv_σmax_upa(mat, state, up)
        dσmaxdΔλ = m*(norm_r + Δλ*dot(r/norm_r, drdΔλ))
        dfdΔλ = dot(dfdσ, dσdΔλ) + dfdσmax*dσmaxdΔλ
        Δλ = Δλ - f/dfdΔλ
        
        if Δλ<=0 || isnan(Δλ) || i==maxits
            # return 0.0, state.σ, 0.0, failure("TCJoint: failed to find Δλ")
            # switch to bissection method
            return calc_σ_up_Δλ_bis(mat, state, σtr)
        end

        if maximum(abs, σ-σ0) <= tol
            break
        end
        σ0 .= σ
    end

    return σ, up, Δλ, success()
end


function calc_σ_up(mat::TCJoint, state::TCJointState, σtr::Array{Float64,1}, Δλ::Float64)
    ndim = state.env.ndim
    kn, ks  = calc_kn_ks(mat, state)

    if ndim == 3
        if σtr[1]>0
            σ = [ σtr[1]/(1+2*Δλ*kn),  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
        else
            σ = [ σtr[1],  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
        end
    else
        if σtr[1]>0
            σ = [ σtr[1]/(1+2*Δλ*kn),  σtr[2]/(1+2*Δλ*ks) ]
        else
            σ = [ σtr[1],  σtr[2]/(1+2*Δλ*ks) ]
        end
    end

    r  = potential_derivs(mat, state, σ)
    up = state.up + Δλ*norm(r)
    return σ, up
end

function calc_σ_up_Δλ_bis(mat::TCJoint, state::TCJointState, σtr::Array{Float64,1})
    ndim    = state.env.ndim
    kn, ks  = calc_kn_ks(mat, state)
    De      = diagm([kn, ks, ks][1:ndim])
    r       = potential_derivs(mat, state, state.σ)

    # Δλ estimative
    Δλ0 = norm(σtr-state.σ)/norm(De*r)
    
    # find initial interval
    # a  = 0.0
    # b  = Δλ0

    ff(Δλ)  = begin
        # quantities at n+1
        σ, up = calc_σ_up(mat, state, σtr, Δλ)
        σmax = calc_σmax(mat, state, up)
        yield_func(mat, state, σ, σmax)
    end

    a, b, status = findrootinterval(ff, 0.0, Δλ0)
    failed(status) && return state.σ, 0.0, 0.0, status
    # @show a, b

    Δλ, status = findroot(ff, a, b, ftol=1e-5, method=:bisection)
    failed(status) && return state.σ, 0.0, 0.0, status

    σ, up = calc_σ_up(mat, state, σtr, Δλ)
    return σ, up, Δλ, success()  

end


function yield_func_from_Δλ(mat::TCJoint, state::TCJointState, σtr::Array{Float64,1}, Δλ::Float64)
    ndim = state.env.ndim
    kn, ks = calc_kn_ks(mat, state)

    # quantities at n+1
    if ndim == 3
        if σtr[1]>0
            σ = [ σtr[1]/(1+2*Δλ*kn),  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
        else
            σ = [ σtr[1],  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
        end
    else
        if σtr[1]>0
            σ = [ σtr[1]/(1+2*Δλ*kn),  σtr[2]/(1+2*Δλ*ks) ]
        else
            σ = [ σtr[1],  σtr[2]/(1+2*Δλ*ks) ]
        end
    end

    r  = potential_derivs(mat, state, σ)
    nr = norm(r)
    up = state.up + Δλ*nr
    
    σmax = calc_σmax(mat, state, up)
    f    = yield_func(mat, state, σ, σmax)

    return f
end


function calc_σ_up_Δλ_bisection(mat::TCJoint, state::TCJointState, σtr::Array{Float64,1})
    ndim    = state.env.ndim
    kn, ks  = calc_kn_ks(mat, state)
    De      = diagm([kn, ks, ks][1:ndim])
    r       = potential_derivs(mat, state, state.σ)

    # Δλ estimative
    Δλ0 = norm(σtr-state.σ)/norm(De*r)
    
    # find initial interval
    a  = 0.0
    b  = Δλ0
    fa = yield_func_from_Δλ(mat, state, σtr, a)
    fb = yield_func_from_Δλ(mat, state, σtr, b)

    # search for a valid interval
    if fa*fb>0
        δ = Δλ0
        maxits = 100
        for i in 1:maxits
            b += δ
            fb = yield_func_from_Δλ(mat, state, σtr, b)
            if fa*fb<0.0
                break
            end
            δ *= 1.5

            if i==maxits
                return 0.0, state.σ, 0.0, failure("TCJoint: could not find iterval for Δλ")
            end
        end
    end

    # bissection method
    Δλ = 0.0
    f  = 0.0
    up = 0.0
    σ0 = zeros(ndim) # previous value
    σ  = zeros(ndim)
    tol = 1e-5
    maxits = 100
    for i in 1:maxits
        Δλ = (a+b)/2

        f = yield_func_from_Δλ(mat, state, σtr, Δλ)
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

        i==maxits && return 0.0, state.σ, 0.0, failure("TCJoint: could not find Δλ with NR/bissection (maxits reached, f=$f)")
    end

    return σ, up, Δλ, success()      
end


function update_state!(mat::TCJoint, state::TCJointState, Δw::Array{Float64,1})

    ndim = state.env.ndim
    σini = copy(state.σ)

    kn, ks = calc_kn_ks(mat, state)
    De = diagm([kn, ks, ks][1:ndim])
    σmax = calc_σmax(mat, state, state.up)  

    if isnan(Δw[1]) || isnan(Δw[2])
        alert("TCJoint: Invalid value for joint displacement: Δw = $Δw")
    end

    # σ trial and F trial
    σtr  = state.σ + De*Δw

    Ftr  = yield_func(mat, state, σtr, σmax)

    # Elastic and EP integration
    # if σmax == 0.0 && state.w[1] >= 0.0
    #     # Return to apex:
    #     if ndim==3
    #         r1 = [ σtr[1]/kn, σtr[2]/ks, σtr[3]/ks ]
    #         r = r1/norm(r1)
    #         state.Δλ = norm(r1)
    #     else
    #         r1 = [ σtr[1]/kn, σtr[2]/ks ]
    #         r = r1/norm(r1)
    #         state.Δλ = norm(r1)  
    #     end

    #     state.up += state.Δλ
    #     state.σ = σtr - state.Δλ*De*r

    # else
    if Ftr <= 0.0
        state.Δλ  = 0.0
        state.σ  .= σtr
    elseif state.up>=mat.wc && σtr[1]>0
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
        σ, up, Δλ, status = calc_σ_up_Δλ(mat, state, σtr)
        failed(status) && return state.σ, status

        state.σ, state.up, state.Δλ = σ, up, Δλ
    end
    state.w += Δw
    Δσ = state.σ - σini
    return Δσ, success()
end


function ip_state_vals(mat::TCJoint, state::TCJointState)
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


function output_keys(mat::TCJoint)
    return Symbol[:jw1, :js1, :jup]
end