 #This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export MMCJoint

mutable struct MMCJointIpState<:IpState
    env::ModelEnv
    σ  ::Array{Float64,1} # stress
    w  ::Array{Float64,1} # relative displacements
    upa::Float64          # effective plastic relative displacement
    Δλ ::Float64          # plastic multiplier
    h  ::Float64          # characteristic length from bulk elements
    function MMCJointIpState(env::ModelEnv=ModelEnv())
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

# Returns the element type that works with this material model
matching_elem_type(::MMCJoint) = MechJoint

# Type of corresponding state structure
ip_state_type(mat::MMCJoint) = MMCJointIpState


function yield_func(mat::MMCJoint, ipd::MMCJointIpState, σ::Array{Float64,1}, σmax::Float64)
    ft, α, β = mat.ft, mat.α, mat.β

    if ipd.env.ndim == 3
        return σ[1] - σmax + β*((σ[2]^2 + σ[3]^2)/ft^2)^α
    else
        return σ[1] - σmax + β*(σ[2]^2/ft^2)^α
    end
end


function yield_derivs(mat::MMCJoint, ipd::MMCJointIpState, σ::Array{Float64,1})
    ft, α, β = mat.ft, mat.α, mat.β

    if ipd.env.ndim == 3
        tmp = 2*α*β/ft^2*((σ[2]^2+σ[3]^2)/ft^2)^(α-1)
        return [ 1 , σ[2]*tmp, σ[3]*tmp ]
    else
        tmp = 2*α*β/ft^2*(σ[2]^2/ft^2)^(α-1)
        return [ 1 , σ[2]*tmp ]
    end
end


function potential_derivs(mat::MMCJoint, ipd::MMCJointIpState, σ::Array{Float64,1})
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


function calc_σmax(mat::MMCJoint, ipd::MMCJointIpState, upa::Float64)
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


function deriv_σmax_upa(mat::MMCJoint, ipd::MMCJointIpState, upa::Float64)
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


function calc_kn_ks(mat::MMCJoint, ipd::MMCJointIpState)
    kn = mat.E*mat.ζ/ipd.h
    G  = mat.E/(2.0*(1.0+mat.ν))
    ks = G*mat.ζ/ipd.h

    return kn, ks
end


function calc_Δλ(mat::MMCJoint, ipd::MMCJointIpState, σtr::Array{Float64,1})
    ndim = ipd.env.ndim
    maxits = 20
    Δλ     = 0.0
    f      = 0.0
    upa    = 0.0
    tol    = 1e-4
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

        f = yield_func(mat, ipd, σ, σmax)
        dfdσ = yield_derivs(mat, ipd, σ)

        m = deriv_σmax_upa(mat, ipd, upa)
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
            return 0.0, failure()
        end
    end
    # @show nits
    # @show f, Δλ

    return Δλ, success()
end


function calc_σ_upa(mat::MMCJoint, ipd::MMCJointIpState, σtr::Array{Float64,1})
    ndim = ipd.env.ndim
    kn, ks = calc_kn_ks(mat, ipd)

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
    return ipd.σ, ipd.upa
end


function mountD(mat::MMCJoint, ipd::MMCJointIpState)
    ndim = ipd.env.ndim
    kn, ks = calc_kn_ks(mat, ipd)
    σmax = calc_σmax(mat, ipd, ipd.upa)

    De = diagm([kn, ks, ks][1:ndim])
    # @show σmax
    # @show ipd.upa
    # @show ipd.upa > mat.wc

    if ipd.Δλ == 0.0  # Elastic 
        # @show "ELASTIC De"
        return De
    # elseif σmax == 0.0 # ?????????????? and w1>0 ???
    elseif σmax == 0.0 && ipd.w[1]>0
        # @show "smax=0 De"

        Dep = De*1e-4
        # Dep  = diagm(ones(ndim))
        # Dep[1,1] *= -1
        return Dep
    else
        # @show "plastic De"
        r = potential_derivs(mat, ipd, ipd.σ)
        v = yield_derivs(mat, ipd, ipd.σ)
        y = -1  # ∂F/∂σmax
        m = deriv_σmax_upa(mat, ipd, ipd.upa)  # ∂σmax/∂upa

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


function stress_update(mat::MMCJoint, ipd::MMCJointIpState, Δw::Array{Float64,1})

    ndim = ipd.env.ndim
    σini = copy(ipd.σ)

    kn, ks = calc_kn_ks(mat, ipd)
    De = diagm([kn, ks, ks][1:ndim])
    σmax = calc_σmax(mat, ipd, ipd.upa)  
    # @show σmax

    if isnan(Δw[1]) || isnan(Δw[2])
        alert("MMCJoint: Invalid value for joint displacement: Δw = $Δw")
    end

    # σ trial and F trial
    σtr  = ipd.σ + De*Δw

    Ftr  = yield_func(mat, ipd, σtr, σmax)
    # @show "stress update"
    # @show σmax
    # @show ipd.upa
    # @show Δw
    # @show ipd.σ
    # @show σtr
    # @show Ftr
    # @show ipd.w[1] 

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

        ipd.σ, ipd.upa = calc_σ_upa(mat, ipd, σtr)

        # @show ipd.Δλ
                      
        # Return to surface:
        # F  = yield_func(mat, ipd, ipd.σ)   
        # F > 1e-2 && alert("MMCJoint: Yield function value ($F) outside tolerance")

    end
    ipd.w += Δw
    Δσ = ipd.σ - σini
    return Δσ, success()
end


function ip_state_vals(mat::MMCJoint, ipd::MMCJointIpState)
    ndim = ipd.env.ndim
    if ndim == 3
       return Dict(
          :wn  => ipd.w[1],
          :w2  => ipd.w[2],
          :w3  => ipd.w[3],
          :sn  => ipd.σ[1],
          :s2  => ipd.σ[2],
          :s3  => ipd.σ[3],
          :upa => ipd.upa
          )
    else
        return Dict(
          :wn  => ipd.w[1],
          :w2  => ipd.w[2],
          :sn  => ipd.σ[1],
          :s2  => ipd.σ[2],
          :upa => ipd.upa
          )
    end
end

function output_keys(mat::MMCJoint)
    return Symbol[:wn, :sn, :upa]
end