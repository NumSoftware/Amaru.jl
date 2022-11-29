 #This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export TCJoint

abstract type AbstractTCJointState<:IpState end

mutable struct TCJointState<:AbstractTCJointState
    env::ModelEnv
    σ  ::Array{Float64,1} # stress
    w  ::Array{Float64,1} # relative displacements
    up ::Float64          # effective plastic relative displacement
    Δλ ::Float64          # plastic multiplier
    h  ::Float64          # characteristic length from bulk elements
    function TCJointState(env::ModelEnv=ModelEnv())
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

abstract type AbstractTCJoint<:Material end

mutable struct TCJoint<:AbstractTCJoint
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

function paramsdict(mat::AbstractTCJoint)
    params = OrderedDict( string(field)=> getfield(mat, field) for field in fieldnames(typeof(mat)) )

    mat.softcurve in ("hordijk", "soft") && ( params["GF"] = 0.1943*mat.ft*mat.wc )
    return params
end


# Returns the element type that works with this material model
matching_elem_type(::AbstractTCJoint) = MechJoint

# Type of corresponding state structure
ip_state_type(::AbstractTCJoint) = TCJointState

function beta(mat::AbstractTCJoint, σmax::Float64)
    βini = mat.βini
    βres = mat.γ*βini
    return βres + (βini-βres)*(σmax/mat.ft)^mat.α
end


function yield_func(mat::AbstractTCJoint, state::AbstractTCJointState, σ::Array{Float64,1}, σmax::Float64)
    α  = mat.α
    β = beta(mat, σmax)
    ft = mat.ft
    if state.env.ndim == 3
        return β*(σ[1] - σmax) + ((σ[2]^2 + σ[3]^2)/ft^2)^α
    else
        return β*(σ[1] - σmax) + (σ[2]^2/ft^2)^α
    end
end


function yield_derivs(mat::AbstractTCJoint, state::AbstractTCJointState, σ::Array{Float64,1}, σmax::Float64)
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


function potential_derivs(mat::AbstractTCJoint, state::AbstractTCJointState, σ::Array{Float64,1})
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


function calc_σmax(mat::AbstractTCJoint, state::AbstractTCJointState, up::Float64)
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
        c = 0.1 
        a = 1.10366
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


function deriv_σmax_upa(mat::AbstractTCJoint, state::AbstractTCJointState, up::Float64)
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
            dz = ((81*up^2*exp(-6.93*up/mat.wc)/mat.wc^3) - (6.93*(1 + 27*up^3/mat.wc^3)*exp(-6.93*up/mat.wc)/mat.wc) - 0.02738402432/mat.wc)
        else
            dz = 0.0
        end
        dσmax = dz*mat.ft 
    elseif mat.softcurve == "soft"
        # c = 1.5
        # a = 1.19311
        # c = 0.0
        # a = 1.09913
        # c = 2.0
        # a = 1.24787
        c = 0.1
        a = 1.10366
        # a = 1.07946
        m = 0.55
        a = 1.30837

        if up == 0.0
            dz = 0.0
        elseif up < mat.wc
            x = up/mat.wc
            # dz = c*(1.0-x)^(c-1)*(a^(1.0-1.0/x)-1.0) - (a^(1.0-1.0/x)*log(a)*(1.0-x)^c)/x^2
            dz =  - (a^(1.0-1.0/x^m)*m*log(a))/x^(m+1)

        else
            dz = 0.0
        end
        dσmax = dz*mat.ft 
    end

    # @show mat.softcurve


    return dσmax
end


function calc_kn_ks(mat::AbstractTCJoint, state::AbstractTCJointState)
    kn = mat.E*mat.ζ/state.h
    G  = mat.E/(2.0*(1.0+mat.ν))
    ks = G*mat.ζ/state.h

    return kn, ks
end


function calc_Δλ(mat::AbstractTCJoint, state::AbstractTCJointState, σtr::Array{Float64,1})
    ndim = state.env.ndim
    maxits = 50
    Δλ     = 0.0
    f      = 0.0
    up    = 0.0
    # tol    = 1e-2 # bad
    # tol    = 1e-4 # better
    tol    = 1e-5 # best
    ft = mat.ft
    θ    = mat.θ
    βini = mat.βini
    βres = mat.γ*βini
    
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
        β      = beta(mat, σmax)

        f    = yield_func(mat, state, σ, σmax)
        dfdσ = yield_derivs(mat, state, σ, σmax)

        dfdσmax = (βini-βres)/ft*(σ[1]-σmax)*θ*(σmax/ft)^(θ-1) - β
        m = deriv_σmax_upa(mat, state, up)
        dσmaxdΔλ = m*(norm_r + Δλ*dot(r/norm_r, drdΔλ))
        dfdΔλ = dot(dfdσ, dσdΔλ) + dfdσmax*dσmaxdΔλ
        Δλ = Δλ - f/dfdΔλ
        
        abs(f) < tol && break

        if i == maxits || isnan(Δλ) || Δλ<=0
            return 0.0, failure("TCJoint: could not find Δλ")
        end
    end

    return Δλ, success()
end


function consistentD(mat::AbstractTCJoint, state::AbstractTCJointState)
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
        Δσ, succeeded = stress_update(mat, statej, Δw)
        # @show Δw
        # @show Δσ
        Dep[:,j] .= Δσ./h
        V[j] = 0.0
    end

    # @show Dep
    # error()

    return Dep
end


function mountD(mat::AbstractTCJoint, state::AbstractTCJointState)
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

        fc, ft = mat.fc, mat.ft
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


function stress_update(mat::AbstractTCJoint, state::AbstractTCJointState, Δw::Array{Float64,1})

    ndim = state.env.ndim
    σini = copy(state.σ)

    kn, ks = calc_kn_ks(mat, state)
    De = diagm([kn, ks, ks][1:ndim])
    σmax = calc_σmax(mat, state, state.up)  

    if isnan(Δw[1]) || isnan(Δw[2])
        alert("AbstractTCJoint: Invalid value for joint displacement: Δw = $Δw")
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
        state.Δλ = 0.0
        state.σ  = copy(σtr) 

    else
        # Plastic increment
        state.Δλ, status = calc_Δλ(mat, state, σtr)
        failed(status) && return state.σ, status

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


function ip_state_vals(mat::AbstractTCJoint, state::AbstractTCJointState)
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

function output_keys(mat::AbstractTCJoint)
    return Symbol[:jw1, :js1, :jup]
end
    