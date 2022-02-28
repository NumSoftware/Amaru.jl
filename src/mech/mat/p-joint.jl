 #This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export PJoint

mutable struct PJointIpState<:IpState
    env::ModelEnv
    σ  ::Array{Float64,1} # stress
    w  ::Array{Float64,1} # relative displacements
    upa::Float64          # effective plastic relative displacement
    Δλ ::Float64          # plastic multiplier
    h  ::Float64          # characteristic length from bulk elements
    function PJointIpState(env::ModelEnv=ModelEnv())
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

# Returns the element type that works with this material model
matching_elem_type(::PJoint) = MechJoint

# Type of corresponding state structure
ip_state_type(mat::PJoint) = PJointIpState

function yield_func(mat::PJoint, ipd::PJointIpState, σ::Array{Float64,1}, σmax::Float64)
    ndim = ipd.env.ndim
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


function yield_derivs(mat::PJoint, ipd::PJointIpState, σ::Array{Float64,1}, σmax::Float64)
    ndim = ipd.env.ndim
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


function potential_derivs(mat::PJoint, ipd::PJointIpState, σ::Array{Float64,1})
    ndim = ipd.env.ndim
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


function calc_σmax(mat::PJoint, ipd::PJointIpState, upa::Float64)
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
            e = exp(1.0)
            z = (1 + 27*(upa/mat.wc)^3)*e^(-6.93*upa/mat.wc) - 28*(upa/mat.wc)*e^(-6.93)
        else
            z = 0.0
        end
        σmax = z*mat.ft 
    end

    # σmax<0.001*mat.ft  && (σmax=0.0)

    return σmax
end


function deriv_σmax_upa(mat::PJoint, ipd::PJointIpState, upa::Float64)
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
            e = exp(1.0)
            dz = ((81*upa^2*e^(-6.93*upa/mat.wc)/mat.wc^3) - (6.93*(1 + 27*upa^3/mat.wc^3)*e^(-6.93*upa/mat.wc)/mat.wc) - 0.02738402432/mat.wc)
        else
            dz = 0.0
        end
        dσmax = dz*mat.ft 
    end

    return dσmax
end


function calc_kn_ks(mat::PJoint, ipd::PJointIpState)
    kn = mat.E*mat.ζ/ipd.h
    G  = mat.E/(2.0*(1.0+mat.ν))
    ks = G*mat.ζ/ipd.h

    return kn, ks
end


function calc_Δλ(mat::PJoint, ipd::PJointIpState, σtr::Array{Float64,1})
    ndim = ipd.env.ndim
    maxits = 20
    Δλ     = 0.0
    f      = 0.0
    upa    = 0.0
    tol    = 1e-2
    fc, ft = mat.fc, mat.ft
    βini = 2*ft - fc -2*√(ft^2 - fc*ft)
    βres = mat.γ*βini
    α    = mat.α
    
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
        m = deriv_σmax_upa(mat, ipd, upa)
        dσmaxdΔλ = m*(norm_r + Δλ*dot(r/norm_r, drdΔλ))
        dfdΔλ = dot(dfdσ, dσdΔλ) + dfdσmax*dσmaxdΔλ
        Δλ = Δλ - f/dfdΔλ
        
        abs(f) < tol && break

        if i == maxits || isnan(Δλ)
            return 0.0, failure()
        end
    end

    return Δλ, success()
end


# function calc_σ_upa(mat::PJoint, ipd::PJointIpState, σtr::Array{Float64,1})
#     ndim = ipd.env.ndim
#     kn, ks = calc_kn_ks(mat, ipd)

#     if ndim == 3
#         if σtr[1] > 0
#             σ = [σtr[1]/(1 + 2*ipd.Δλ*kn), σtr[2]/(1 + 2*ipd.Δλ*ks), σtr[3]/(1 + 2*ipd.Δλ*ks)]
#         else
#             σ = [σtr[1], σtr[2]/(1 + 2*ipd.Δλ*ks), σtr[3]/(1 + 2*ipd.Δλ*ks)]
#         end    
#     else
#         if σtr[1] > 0
#             σ = [σtr[1]/(1 + 2*ipd.Δλ*kn), σtr[2]/(1 + 2*ipd.Δλ*ks)]
#         else
#             σ = [σtr[1], σtr[2]/(1 + 2*ipd.Δλ*ks)]
#         end    
#     end
#     r = potential_derivs(mat, ipd, σ)
#     upa = ipd.upa + ipd.Δλ*norm(r)
#     return σ, upa
# end


function mountD(mat::PJoint, ipd::PJointIpState)
    ndim = ipd.env.ndim
    kn, ks = calc_kn_ks(mat, ipd)
    α = mat.α
    σmax = calc_σmax(mat, ipd, ipd.upa)

    De = diagm([kn, ks, ks][1:ndim])
    # @show σmax
    # @show ipd.upa
    # @show ipd.upa > mat.wc

    if ipd.Δλ == 0.0  # Elastic 
        # @show "ELASTIC De"
        return De
    elseif σmax == 0.0 && ipd.w[1] >= 0.0
        # @show "smax=0 De"

        Dep = De*1e-4
        # Dep  = diagm(ones(ndim))
        # Dep[1,1] *= -1
        return Dep
    else
        # @show "plastic De"
        fc, ft = mat.fc, mat.ft
        βini = 2*ft - fc -2*√(ft^2 - fc*ft)
        βres = mat.γ*βini

        r = potential_derivs(mat, ipd, ipd.σ)
        v = yield_derivs(mat, ipd, ipd.σ, σmax)
        # dfdσmax = (βres-βini)/ft*(2*σmax-ipd.σ[1]) - βres
        β = βres + (βini-βres)*(σmax/ft)^mat.α
        dfdσmax = (βini-βres)/ft*(ipd.σ[1]-σmax)*α*(σmax/ft)^(α-1) - β
        # dfdσmax = -β  # ∂F/∂σmax
        m = deriv_σmax_upa(mat, ipd, ipd.upa)  # ∂σmax/∂upa

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


        # if σmax == 0.0 
        #     # @show Dep*den
        #     # Dep += 1e-8*De
        #     # @show De
        #     @show kn
        #     @show ks
        #     @show r
        #     @show v
        #     @show dfdσmax
        #     @show m
        #     @show den
        #     @show Dep
        #     error()
        # end

        return Dep
    end
end


function stress_update(mat::PJoint, ipd::PJointIpState, Δw::Array{Float64,1})

    ndim = ipd.env.ndim
    σini = copy(ipd.σ)

    kn, ks = calc_kn_ks(mat, ipd)
    De = diagm([kn, ks, ks][1:ndim])
    σmax = calc_σmax(mat, ipd, ipd.upa)  
    # @show σmax

    if isnan(Δw[1]) || isnan(Δw[2])
        alert("PJoint: Invalid value for joint displacement: Δw = $Δw")
    end

    # σ trial and F trial
    σtr  = ipd.σ + De*Δw

    Ftr  = yield_func(mat, ipd, σtr, σmax)

    # Elastic and EP integration
    if σmax == 0.0 && ipd.w[1] >= 0.0
        # @show "smax=0 up"
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
        # @show "ELASTIC up"
        # Pure elastic increment
        ipd.Δλ = 0.0
        ipd.σ  = copy(σtr) 

    else
        # @show "plastic up"

        # Plastic increment
        ipd.Δλ, status = calc_Δλ(mat, ipd, σtr)
        failed(status) && return ipd.σ, status

        # ipd.σ, ipd.upa = calc_σ_upa(mat, ipd, σtr)
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


function ip_state_vals(mat::PJoint, ipd::PJointIpState)
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

function output_keys(mat::PJoint)
    return Symbol[:jw1, :js1, :jup]
end