 #This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export TCJoint

mutable struct TCJointIpState<:IpState
    env::ModelEnv
    σ  ::Array{Float64,1} # stress
    w  ::Array{Float64,1} # relative displacements
    upa::Float64          # effective plastic relative displacement
    Δλ ::Float64          # plastic multiplier
    h  ::Float64          # characteristic length from bulk elements
    function TCJointIpState(env::ModelEnv=ModelEnv())
        this = new(env)
        ndim = env.ndim
        this.σ = zeros(ndim)
        this.w = zeros(ndim)
        this.upa = 0.0
        this.Δλ  = 0.0
        this.h  = 0.0
        return this
    end
end

mutable struct TCJoint<:Material
    E ::Float64      # Young's modulus
    ν ::Float64      # Poisson ratio
    ft::Float64      # tensile strength (internal variable)
    fc::Float64      # compressive strength (internal variable)
    ζ ::Float64      # factor ζ controls the elastic relative displacements (formerly α)
    wc::Float64      # critical crack opening
    ws::Float64      # openning at inflection (where the curve slope changes)
    softcurve::String # softening curve model ("linear" or bilinear" or "hordijk")
    α::Float64       # curvature coefficient
    γ::Float64       # factor for βres
    θ::Float64       # fator for the surface reduction speed
    βini::Float64    # initial curvature size

    function TCJoint(prms::Dict{Symbol,Float64})
        return  TCJoint(;prms...)
    end

    function TCJoint(;E=NaN, nu=NaN, fc=NaN, ft=NaN, zeta=NaN, alpha=1.0, gamma=1.0, theta=1.0, wc=NaN, ws=NaN, GF=NaN, Gf=NaN, softcurve="bilinear")

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

        @check E>0.0    
        @check 0<=nu<0.5
        @check fc<0  
        @check ft>=0  
        @check zeta>0
        @check wc>0  
        @check isnan(ws) || ws>0
        softcurve in ("linear", "bilinear", "hordijk") || error("Invalid softcurve: softcurve must to be linear or bilinear or hordijk")

        a = (2*alpha*ft + alpha*fc - fc - √(alpha^2*fc^2 - 4*alpha^2*fc*ft + 4*alpha^2*ft^2 - 2*alpha*fc^2 + fc^2)) / (4*alpha-2)
        b = √(alpha*(2*a-fc)*(ft-a))
        βini = (b^2/ft^2)^alpha/(ft-a)

        this = new(E, nu, ft, fc, zeta, wc, ws, softcurve, alpha, gamma, theta, βini)
        return this
    end
end

# Returns the element type that works with this material model
matching_elem_type(::TCJoint) = MechJoint

# Type of corresponding state structure
ip_state_type(mat::TCJoint) = TCJointIpState

function beta(mat::TCJoint, σmax::Float64)
    βini = mat.βini
    βres = mat.γ*βini
    return βres + (βini-βres)*(σmax/mat.ft)^mat.α
end


function yield_func(mat::TCJoint, ipd::TCJointIpState, σ::Array{Float64,1}, σmax::Float64)
    α  = mat.α
    β = beta(mat, σmax)
    ft = mat.ft
    if ipd.env.ndim == 3
        return β*(σ[1] - σmax) + ((σ[2]^2 + σ[3]^2)/ft^2)^α
    else
        return β*(σ[1] - σmax) + (σ[2]^2/ft^2)^α
    end
end


function yield_derivs(mat::TCJoint, ipd::TCJointIpState, σ::Array{Float64,1}, σmax::Float64)
    α = mat.α
    β = beta(mat, σmax)
    ft = mat.ft

    # tmp = 2*α/ft^2*(σ[2]^2/ft^2)^(α-1)
    
    # isnan(σ[2]) && @show σ[2]
    # isnan(tmp) && @show tmp
    # if isnan(σ[2]*tmp) 
    #     @show tmp
    #     @show σ
    #     @show α
    #     @show ft
    #     @show ft
    # end

    if ipd.env.ndim == 3
        tmp = 2*α/ft^2*((σ[2]^2+σ[3]^2)/ft^2)^(α-1)
        σ[2]==σ[3]==0.0 && (tmp=0)
        return [ β , σ[2]*tmp, σ[3]*tmp ]
    else
        tmp = 2*α/ft^2*(σ[2]^2/ft^2)^(α-1)
        σ[2]==0.0 && (tmp=0)
        return [ β , σ[2]*tmp ]
    end
end


function potential_derivs(mat::TCJoint, ipd::TCJointIpState, σ::Array{Float64,1})
    ndim = ipd.env.ndim
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


function calc_σmax(mat::TCJoint, ipd::TCJointIpState, upa::Float64)
    if mat.softcurve == "linear"
        if upa < mat.wc
            a = mat.ft 
            b = mat.ft /mat.wc
        else
            a = 0.0
            b = 0.0
        end
        σmax = a - b*upa
    elseif mat.softcurve == "bilinear"
        σs = 0.25*mat.ft 
        if upa < mat.ws
            a  = mat.ft  
            b  = (mat.ft  - σs)/mat.ws
        elseif upa < mat.wc
            a  = mat.wc*σs/(mat.wc-mat.ws)
            b  = σs/(mat.wc-mat.ws)
        else
            a = 0.0
            b = 0.0
        end
        σmax = a - b*upa
    elseif mat.softcurve == "hordijk"
        if upa < mat.wc
            z = (1 + 27*(upa/mat.wc)^3)*exp(-6.93*upa/mat.wc) - 28*(upa/mat.wc)*exp(-6.93)
        else
            z = 0.0
        end
        σmax = z*mat.ft 
    end

    return σmax
end


function deriv_σmax_upa(mat::TCJoint, ipd::TCJointIpState, upa::Float64)
    # ∂σmax/∂upa = dσmax
    if mat.softcurve == "linear"
        if upa < mat.wc
            b = mat.ft /mat.wc
        else
            b = 0.0
        end
        dσmax = -b
    elseif mat.softcurve == "bilinear"
        σs = 0.25*mat.ft 
        if upa < mat.ws
            b  = (mat.ft  - σs)/mat.ws
        elseif upa < mat.wc
            b  = σs/(mat.wc-mat.ws)
        else
            b = 0.0
        end
        dσmax = -b
    elseif mat.softcurve == "hordijk"
        if upa < mat.wc
            dz = ((81*upa^2*exp(-6.93*upa/mat.wc)/mat.wc^3) - (6.93*(1 + 27*upa^3/mat.wc^3)*exp(-6.93*upa/mat.wc)/mat.wc) - 0.02738402432/mat.wc)
        else
            dz = 0.0
        end
        dσmax = dz*mat.ft 
    end

    return dσmax
end


function calc_kn_ks(mat::TCJoint, ipd::TCJointIpState)
    kn = mat.E*mat.ζ/ipd.h
    G  = mat.E/(2.0*(1.0+mat.ν))
    ks = G*mat.ζ/ipd.h

    return kn, ks
end


function calc_Δλ(mat::TCJoint, ipd::TCJointIpState, σtr::Array{Float64,1})
    ndim = ipd.env.ndim
    maxits = 30
    Δλ     = 0.0
    f      = 0.0
    upa    = 0.0
    tol    = 1e-4
    tol    = 1e-3
    tol    = 1e-2
    ft = mat.ft
    θ    = mat.θ
    βini = mat.βini
    βres = mat.γ*βini
    
    nits = 0

    for i in 1:maxits
        nits = i
        kn, ks = calc_kn_ks(mat, ipd)

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
                 
        r      = potential_derivs(mat, ipd, σ)
        norm_r = norm(r)
        upa    = ipd.upa + Δλ*norm_r
        σmax   = calc_σmax(mat, ipd, upa)
        β      = beta(mat, σmax)

        f    = yield_func(mat, ipd, σ, σmax)
        dfdσ = yield_derivs(mat, ipd, σ, σmax)

        dfdσmax = (βini-βres)/ft*(σ[1]-σmax)*θ*(σmax/ft)^(θ-1) - β
        m = deriv_σmax_upa(mat, ipd, upa)
        dσmaxdΔλ = m*(norm_r + Δλ*dot(r/norm_r, drdΔλ))
        dfdΔλ = dot(dfdσ, dσdΔλ) + dfdσmax*dσmaxdΔλ
        Δλ = Δλ - f/dfdΔλ
        
        abs(f) < tol && break

        if i == maxits || isnan(Δλ)
            @show i, Δλ
            return 0.0, failure()
        end
    end

    return Δλ, success()
end


function mountD(mat::TCJoint, ipd::TCJointIpState)
    ndim = ipd.env.ndim
    kn, ks = calc_kn_ks(mat, ipd)
    θ = mat.θ
    σmax = calc_σmax(mat, ipd, ipd.upa)

    De = diagm([kn, ks, ks][1:ndim])

    if ipd.Δλ == 0.0  # Elastic 
        return De
    elseif σmax == 0.0 && ipd.w[1] >= 0.0
        Dep = De*1e-4
        # Dep = De*1e-3
        # Dep = De*1e-2
        # Dep = De*1e-2
        # Dep = De*1e-1
        # Dep = De
        return Dep
    else
        fc, ft = mat.fc, mat.ft
        βini = mat.βini
        βres = mat.γ*βini
        β = beta(mat, σmax)
        dfdσmax = (βini-βres)/ft*(ipd.σ[1]-σmax)*θ*(σmax/ft)^(θ-1) - β

        r = potential_derivs(mat, ipd, ipd.σ)
        v = yield_derivs(mat, ipd, ipd.σ, σmax)
        m = deriv_σmax_upa(mat, ipd, ipd.upa)  # ∂σmax/∂upa

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

        if any(isnan.(Dep))
            @show den
            @show ipd.σ
            @show σmax
            @show dfdσmax
            @show m
            @show r[1]
            @show r[2]
            @show v
            @show Dep
        end

        return Dep
    end
end


function stress_update(mat::TCJoint, ipd::TCJointIpState, Δw::Array{Float64,1})

    ndim = ipd.env.ndim
    σini = copy(ipd.σ)

    kn, ks = calc_kn_ks(mat, ipd)
    De = diagm([kn, ks, ks][1:ndim])
    σmax = calc_σmax(mat, ipd, ipd.upa)  

    if isnan(Δw[1]) || isnan(Δw[2])
        alert("TCJoint: Invalid value for joint displacement: Δw = $Δw")
    end

    # σ trial and F trial
    σtr  = ipd.σ + De*Δw

    Ftr  = yield_func(mat, ipd, σtr, σmax)

    # Elastic and EP integration
    if σmax == 0.0 && ipd.w[1] >= 0.0
        # Return to apex:
        if ndim==3
            r1 = [ σtr[1]/kn, σtr[2]/ks, σtr[3]/ks ]
            r = r1/norm(r1)
            ipd.Δλ = norm(r1)
        else
            r1 = [ σtr[1]/kn, σtr[2]/ks ]
            r = r1/norm(r1)
            ipd.Δλ = norm(r1)  
        end

        ipd.upa += ipd.Δλ
        ipd.σ = σtr - ipd.Δλ*De*r

    elseif Ftr <= 0.0
        ipd.Δλ = 0.0
        ipd.σ  = copy(σtr) 

    else
        # Plastic increment
        ipd.Δλ, status = calc_Δλ(mat, ipd, σtr)
        failed(status) && return ipd.σ, status

        if ndim == 3
            if σtr[1] > 0
                ipd.σ = [σtr[1]/(1 + 2*ipd.Δλ*kn), σtr[2]/(1 + 2*ipd.Δλ*ks), σtr[3]/(1 + 2*ipd.Δλ*ks)]
            else
                ipd.σ = [σtr[1], σtr[2]/(1 + 2*ipd.Δλ*ks), σtr[3]/(1 + 2*ipd.Δλ*ks)]
            end    
        else
            if σtr[1] > 0
                ipd.σ = [σtr[1]/(1 + 2*ipd.Δλ*kn), σtr[2]/(1 + 2*ipd.Δλ*ks)]
            else
                ipd.σ = [σtr[1], σtr[2]/(1 + 2*ipd.Δλ*ks)]
            end    
        end
        r = potential_derivs(mat, ipd, ipd.σ)
        ipd.upa += ipd.Δλ*norm(r)

    end
    ipd.w += Δw
    Δσ = ipd.σ - σini
    return Δσ, success()
end


function ip_state_vals(mat::TCJoint, ipd::TCJointIpState)
    ndim = ipd.env.ndim
    if ndim == 3
       return Dict(
          :jw1  => ipd.w[1],
          :jw2  => ipd.w[2],
          :jw3  => ipd.w[3],
          :js1  => ipd.σ[1],
          :js2  => ipd.σ[2],
          :js3  => ipd.σ[3],
          :jup => ipd.upa
          )
    else
        return Dict(
          :jw1  => ipd.w[1],
          :jw2  => ipd.w[2],
          :js1  => ipd.σ[1],
          :js2  => ipd.σ[2],
          :jup => ipd.upa
          )
    end
end

function output_keys(mat::TCJoint)
    return Symbol[:jw1, :js1, :jup]
end
    