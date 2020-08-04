# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export SmearedCrack

mutable struct SmearedCrackIpState<:IpState
    env::ModelEnv
    σ::Tensor2
    ε::Tensor2
    T::Tensor4    # rotation tensor to the crack plane
    w::Array{Float64,1} # relative displacements in the crack plane
    upa::Float64  # max plastic displacement
    Δλ ::Float64  # plastic multiplier
    h  ::Float64  # element size fraction for a integration point
    hascrack::Bool # flag to signal if an ip has cracked
    function SmearedCrackIpState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ = zeros(6)
        this.ε = zeros(6)
        this.T = eye(6)
        this.w   = zeros(3)
        this.upa = 0.0
        this.Δλ  = 0.0
        this.h   = 0.0 # length related with the element size
        this.hascrack = false
        return this
    end
end

mutable struct SmearedCrack<:Material
    E  ::Float64  # Young's modulus from bulk material
    ν  ::Float64  # Poisson ratio from bulk material
    ft::Float64  # tensile strength
    μ  ::Float64  # friction angle
    ζ  ::Float64  # elastic displacement scale factor
    wc ::Float64  # critical openning
    ws ::Float64  # openning at inflection
    softcurve::String # softening model ("bilinear" or "hordijk")

    function SmearedCrack(prms::Dict{Symbol,Float64})
        return SmearedCrack(;prms...)
    end

    function SmearedCrack(;E=NaN, nu=NaN, ft=NaN, mu=NaN, zeta=NaN, wc=NaN, ws=NaN, GF=NaN, Gf=NaN, softcurve="bilinear")
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

        E>0.0       || error("Invalid value for E: $E")
        0<=nu<0.5   || error("Invalid value for nu: $nu")
        ft>=0       || error("Invalid value for ft: $ft")
        mu>0     || error("Invalid value for mu: $mu")
        zeta>0      || error("Invalid value for zeta: $zeta")
        wc>0        || error("Invalid value for wc: $wc")
        (isnan(ws)  || ws>0) || error("Invalid value for ws: $ws")
        (softcurve=="linear" || softcurve=="bilinear" || softcurve=="hordijk") || error("Invalid softcurve: softcurve must to be linear or bilinear or hordijk")

        this = new(E, nu, ft, mu, zeta, wc, ws, softcurve)
        return this
    end
end

# Returns the element type that works with this material model
matching_elem_type(::SmearedCrack) = MechSolid

# Type of corresponding state structure
ip_state_type(mat::SmearedCrack) = SmearedCrackIpState


function calc_σmax(mat::SmearedCrack, ipd::SmearedCrackIpState, upa::Float64)
    wc  = mat.wc
    ws  = mat.ws
    if mat.softcurve == "linear"
        if upa < wc
            a = mat.ft
            b = mat.ft/wc
        else
            a = b = 0.0
        end
        σmax = a - b*upa
    elseif mat.softcurve == "bilinear"
        σs = 0.25*mat.ft
        if upa < ws
            a  = mat.ft
            b  = (mat.ft - σs)/ws
        elseif upa < wc
            a  = wc*σs/(wc-ws)
            b  = σs/(wc-ws)
        else
            a = b = 0.0
        end
        σmax = a - b*upa
    elseif mat.softcurve == "hordijk"
        if upa < wc
            z = (1 + 27*(upa/wc)^3)*exp(-6.93*upa/wc) - 28*(upa/wc)*exp(-6.93)
        else
            z = 0.0
        end
        σmax = z*mat.ft
    end
    return σmax
end


function σmax_deriv(mat::SmearedCrack, ipd::SmearedCrackIpState, upa::Float64)
    # ∂σmax/∂upa = dσmax
    wc  = mat.wc
    ws  = mat.ws
    if mat.softcurve == "linear"
        if upa < wc
            b = mat.ft/wc
        else
            b = 0.0
        end
        dσmax = -b
    elseif mat.softcurve == "bilinear"
        σs = 0.25*mat.ft
        if upa < ws
            b  = (mat.ft - σs)/ws
        elseif upa < wc
            b  = σs/(wc-ws)
        else
            b = 0.0
        end
        dσmax = -b
    elseif mat.softcurve == "hordijk"
        e = exp(1.0)
        if upa < wc
            dz = ((81*upa^2*e^(-6.93*upa/wc)/wc^3) - (6.93*(1 + 27*upa^3/wc^3)*e^(-6.93*upa/wc)/wc) - 0.02738402432/wc)
        else
            dz = 0.0
        end
        dσmax = dz*mat.ft
    end
    return dσmax
end


function yield_func(mat::SmearedCrack, ipd::SmearedCrackIpState, ς::Array{Float64,1})
    σmax = calc_σmax(mat, ipd, ipd.upa)
    return sqrt(ς[2]^2 + ς[3]^2) + (ς[1]-σmax)*mat.μ
end


function yield_deriv(mat::SmearedCrack, ipd::SmearedCrackIpState, ς::Array{Float64,1})
    σ = ς
    τ = sqrt(ipd.σ[2]^2 + ipd.σ[3]^2)
    τ==0.0 && return [ mat.μ, 0.0, 0.0 ]
    return [ mat.μ, ipd.σ[2]/sqrt(ipd.σ[2]^2 + ipd.σ[3]^2), ipd.σ[3]/sqrt(ipd.σ[2]^2 + ipd.σ[3]^2)]
end


function potential_derivs(mat::SmearedCrack, ipd::SmearedCrackIpState, ς::Array{Float64,1})
    σ = ς
    if σ[1] >= 0.0 
        # G1:
        r = [ 2.0*σ[1]*mat.μ^2, 2.0*σ[2], 2.0*σ[3]]
    else
        # G2:
        r = [ 0.0, 2.0*σ[2], 2.0*σ[3] ]
    end
    return r
end


function calc_kn_ks(mat::SmearedCrack, ipd::SmearedCrackIpState)
    kn = mat.E*mat.ζ/ipd.h
    G  = mat.E/(2.0*(1.0+mat.ν))
    ks = G*mat.ζ/ipd.h

    return kn, ks
end


function calc_Δλ(mat::SmearedCrack, ipd::SmearedCrackIpState, σtr::Array{Float64,1})
    maxits = 100
    Δλ     = 0.0
    f      = 0.0
    upa    = 0.0
    tol    = 1e-4      

    for i=1:maxits
        μ      = mat.μ
        kn, ks = calc_kn_ks(mat, ipd)

        # quantities at n+1
        if σtr[1]>0
             σ     = [ σtr[1]/(1+2*Δλ*kn*μ^2),  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
             dσdΔλ = [ -2*kn*μ^2*σtr[1]/(1+2*Δλ*kn*μ^2)^2,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2,  -2*ks*σtr[3]/(1+2*Δλ*ks)^2 ]
             drdΔλ = [ -4*kn*μ^4*σtr[1]/(1+2*Δλ*kn*μ^2)^2,  -4*ks*σtr[2]/(1+2*Δλ*ks)^2,  -4*ks*σtr[3]/(1+2*Δλ*ks)^2 ]
        else
             σ     = [ σtr[1],  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
             dσdΔλ = [ 0,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2,  -2*ks*σtr[3]/(1+2*Δλ*ks)^2 ]
             drdΔλ = [ 0,  -4*ks*σtr[2]/(1+2*Δλ*ks)^2,  -4*ks*σtr[3]/(1+2*Δλ*ks)^2 ]
        end

        r      = potential_derivs(mat, ipd, σ)
        norm_r = norm(r)
        upa    = ipd.upa + Δλ*norm_r
        σmax   = calc_σmax(mat, ipd, upa)
        m      = σmax_deriv(mat, ipd, upa)
        dσmaxdΔλ = m*(norm_r + Δλ*dot(r/norm_r, drdΔλ))

        f = sqrt(σ[2]^2 + σ[3]^2) + (σ[1]-σmax)*μ
        if (σ[2]==0 && σ[3]==0) 
            dfdΔλ = (dσdΔλ[1] - dσmaxdΔλ)*μ
        else
            dfdΔλ = 1/sqrt(σ[2]^2 + σ[3]^2) * (σ[2]*dσdΔλ[2] + σ[3]*dσdΔλ[3]) + (dσdΔλ[1] - dσmaxdΔλ)*μ
        end

        Δλ = Δλ - f/dfdΔλ

        abs(f) < tol && break

        if i == maxits || isnan(Δλ)
            @error """MCJoint: Could not find Δλ. This may happen when the system
            becomes hypostatic and thus the global stiffness matrix is near syngular.
            Increasing the mesh refinement may result in a nonsingular matrix.
            """ iterations=i Δλ
            error()
        end
    end
    return Δλ
end


function calc_σ_upa(mat::SmearedCrack, ipd::SmearedCrackIpState, σtr::Array{Float64,1})
    μ = mat.μ
    kn, ks = calc_kn_ks(mat, ipd)

    if σtr[1] > 0
        σ = [σtr[1]/(1 + 2*ipd.Δλ*kn*(μ^2)), σtr[2]/(1 + 2*ipd.Δλ*ks), σtr[3]/(1 + 2*ipd.Δλ*ks)]
    else
        σ = [σtr[1], σtr[2]/(1 + 2*ipd.Δλ*ks), σtr[3]/(1 + 2*ipd.Δλ*ks)]
    end    
    r = potential_derivs(mat, ipd, σ)
    ipd.upa += ipd.Δλ*norm(r)
    return σ, ipd.upa
end


function calcD(mat::SmearedCrack, ipd::SmearedCrackIpState)
    E, ν = mat.E, mat.ν
    D = calcDe(E, ν, ipd.env.modeltype)

    # If has no crack return complete elastic tensor
    ipd.hascrack || return D

    σmax = calc_σmax(mat, ipd, ipd.upa)
    kn, ks = calc_kn_ks(mat, ipd)

    if ipd.Δλ == 0.0 # Cracked but in elastic regime
        #@show "elastic"
        #if ipd.upa>mat.wc
            #error("here")
        #end
        # Modify terms in the elastic tensor
        #z = kn*ipd.h/mat.ζ*1e-6
        #D[2,1] = D[3,1] = z  #!!!!!!!!!!!!!!!!
        #D[2,3] = D[3,2] = z

        Dcr =  ipd.h/mat.ζ*[  kn 0  0
                              0  ks 0 
                              0  0  ks ]
        Dcr = inv(Dcr)

            D[1,2] = D[1,3] = 0.0
            D[2,1] = D[3,1] = 0.0
            D[2,3] = D[3,2] = 0.0  # sometimes not conv if not zero
            
        #D = inv(D)
        #D[1,1] = Dcr[1,1]
        #D[5,5] = Dcr[2,2]*2
        #D[6,6] = Dcr[3,3]*2

        #return ipd.T'*inv(D)*ipd.T

        D[1,:] .= ( kn*ipd.h/mat.ζ, 0.0, 0.0, 0.0, 0.0, 0.0 )
        D[5,:] .= (        0.0, 0.0, 0.0, 0.0, 2*ks*ipd.h/mat.ζ, 0.0 )
        D[6,:] .= (        0.0, 0.0, 0.0, 0.0, 0.0, 2*ks*ipd.h/mat.ζ )

    else
        if σmax == 0.0 && ipd.w[1] >= 0.0
            #@show "fully cracked"
            #error()
            k = kn*ipd.h/mat.ζ*1e-6
            z = kn*ipd.h/mat.ζ*1e-6
                D[1,2] = D[1,3] = 0.0
                D[2,1] = D[3,1] = 0.0
                D[2,3] = D[3,2] = 0.0

            #display(D)

            D = inv(D)
            k=1\k
            D[1,:] .+= (   k, 0.0, 0.0, 0.0, 0.0, 0.0 )
            D[5,:] .+= ( 0.0, 0.0, 0.0, 0.0,   k, 0.0 )
            D[6,:] .+= ( 0.0, 0.0, 0.0, 0.0, 0.0,   k )
            #display(inv(D))
            #error()
            return ipd.T'*inv(D)*ipd.T

        else # Elastoplastic tensor D
            #@show "partially cracked"
            idx  = [1, 5, 6]          # indexes for components in the crack plane
            mand = (1.0, SR2, SR2)    # Mandel's notation correction coefficients
            ς = (ipd.T*ipd.σ)[idx]./mand  # stress components in the crack plane

            v = yield_deriv(mat, ipd, ς)
            r = potential_derivs(mat, ipd, ς)
            y = -mat.μ # ∂F/∂σmax
            m = σmax_deriv(mat, ipd, ipd.upa)  # ∂σmax/∂upa
            #@show v

            den = kn*r[1]*v[1] + ks*r[2]*v[2] + ks*r[3]*v[3] - y*m*norm(r)
            Dcr =  ipd.h*[  kn - kn^2*r[1]*v[1]/den    -kn*ks*r[1]*v[2]/den      -kn*ks*r[1]*v[3]/den
                           -kn*ks*r[2]*v[1]/den         ks - ks^2*r[2]*v[2]/den  -ks^2*r[2]*v[3]/den
                           -kn*ks*r[3]*v[1]/den        -ks^2*r[3]*v[2]/den        ks - ks^2*r[3]*v[3]/den ]

            # Modify terms in the elastic tensor considering Mandel's notation
            #display(Dcr)
            #display(D)
            z = kn*ipd.h/mat.ζ*1e-6
                D[1,2] = D[1,3] = 0.0
                D[2,1] = D[3,1] = 0.0
                D[2,3] = D[3,2] = 0.0
            D = inv(D)
            #@show Dcr
            #@show den
            Dcr = pinv(Dcr)

            D[1,:] .+= ( Dcr[1,1]    , 0.0, 0.0, 0.0, Dcr[1,2]*SR2, Dcr[1,3]*SR2 )
            D[5,:] .+= ( Dcr[2,1]*SR2, 0.0, 0.0, 0.0, Dcr[2,2]*2    , Dcr[2,3]*2 )
            D[6,:] .+= ( Dcr[3,1]*SR2, 0.0, 0.0, 0.0, Dcr[3,2]*2    , Dcr[3,3]*2 )
            #display(D)
            #display(inv(D))
            #error()

            return ipd.T'*inv(D)*ipd.T

        end
    end

    # Rotate modified D to the global xyz system
    D = ipd.T'*D*ipd.T
    return D
end


function stress_update(mat::SmearedCrack, ipd::SmearedCrackIpState, Δε::Array{Float64,1})
    # σ : stress tensor in the xyz system
    # ε : strain tensor in the xyz system
    # ς : stress vector at crack plane

    E, ν = mat.E, mat.ν
    T = ipd.T
    h = ipd.h
    kn, ks = calc_kn_ks(mat, ipd)
    mand = (1.0, SR2, SR2)
    idx  = [1, 5, 6]     # indexes for components in the crack plane
    D = calcD(mat, ipd)

    σini  = ipd.σ
    # Check for cracks
    if !ipd.hascrack
        Δσtr  = D*Δε
        σtr   = ipd.σ + Δσtr
        P = eigvals(σtr)
        if P[1]>mat.ft
            P, V = eigen(σtr)
            ipd.T = tensor_rot(V)
            ipd.hascrack = true
        else
            ipd.σ = σtr
            return Δσtr
        end
    end

    Decr = [kn/mat.ζ, ks/mat.ζ, ks/mat.ζ]
    ςini = (T*σini)[idx]./mand # initial stress components in the crack plane
    Δw   = (T*Δε)[idx]./mand*ipd.h
    ςtr  = ςini .+ Decr.*Δw
    Ftr  = yield_func(mat, ipd, ςtr)

    σmax = calc_σmax(mat, ipd, ipd.upa)
    if σmax == 0.0 && ipd.w[1] >= 0.0
        r1 = [ ςtr[1]/kn, ςtr[2]/ks, ςtr[3]/ks ]
        ipd.Δλ   = norm(r1)
        r = normalize!(r1)
        ipd.upa += ipd.Δλ
        ς = ςtr - ipd.Δλ*Decr.*r
    elseif Ftr <= 0
        #if ipd.upa > 0 
            #println()
            #Fini  = yield_func(mat, ipd, ςini)
            #@show Fini
            #@show ςtr
            #@show Δw
            #@show Δε
        #end
        ipd.Δλ = 0.0
        ς = ςtr
    else
        ipd.Δλ = calc_Δλ(mat, ipd, ςtr) 
        ς, ipd.upa = calc_σ_upa(mat, ipd, ςtr)

        # Return to surface:
        F  = yield_func(mat, ipd, ς)   
        F > 1e-3 && @warn "SmearedCrack: Yield function value outside tolerance:" F
    end

    Δς = ς - ςini
    Δεcr = T*Δε
    Δσcr = D*Δεcr


    # Fix Δσcr using crack stresses
    Δσcr[idx] .= Δς.*mand
    #@show Δς
    #@show Δσcr

    # Update state
    ipd.ε += Δε
    ipd.w += Δw

    ipd.σ  = σini + T'*Δσcr

    #if σmax == 0.0 && ipd.w[1] >= 0.0
        #σcr = T*σini + Δσcr
        #σcr[idx] .= 0.0
        #ipd.σ  = T'*σcr
    #end

    Δσ     = ipd.σ - σini

    return Δσ
end

function ip_state_vals(mat::SmearedCrack, ipd::SmearedCrackIpState)
    ndim  = ipd.env.ndim
    σ, ε  = ipd.σ, ipd.ε

    D = stress_strain_dict(σ, ε, ipd.env)
    mand = (1.0, SR2, SR2)
    T    = ipd.T
    if ipd.hascrack
        s1 = (T*ipd.σ)[1]
        s2 = (T*ipd.σ)[5]
        s3 = (T*ipd.σ)[6]
    else
        s1 = s2 = s3 = 0.0
    end
    D[:upa] = ipd.upa
    D[:s1]  = s1
    D[:s2]  = s2
    D[:s3]  = s3
    D[:w]   = ipd.w[1]
    D[:w2]  = ipd.w[2]
    D[:w3]  = ipd.w[3]
    D[:crack]  = float(ipd.hascrack)
    return D
end
