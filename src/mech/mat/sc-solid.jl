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
    function SmearedCrackIpState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ = zeros(6)
        this.ε = zeros(6)
        this.T = eye(6)
        this.w   = zeros(3)
        this.upa = 0.0
        this.Δλ  = 0.0
        this.h   = 0.0 # length related with the element size
        return this
    end
end

mutable struct SmearedCrack<:Material
    E  ::Float64  # Young's modulus from bulk material
    ν  ::Float64  # Poisson ratio from bulk material
    ft ::Float64  # tensile strength
    μ  ::Float64  # friction angle
    wc ::Float64  # critical openning
    ws ::Float64  # openning at inflection
    softcurve::String # softening model ("bilinear" or "hordijk")

    function SmearedCrack(prms::Dict{Symbol,Float64})
        return SmearedCrack(;prms...)
    end

    function SmearedCrack(;E=NaN, nu=NaN, ft=NaN, mu=NaN, wc=NaN, ws=NaN, GF=NaN, Gf=NaN, softcurve="hordijk")
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
        wc>0        || error("Invalid value for wc: $wc")
        (isnan(ws)  || ws>0) || error("Invalid value for ws: $ws")
        (softcurve=="linear" || softcurve=="bilinear" || softcurve=="hordijk") || error("Invalid softcurve: softcurve must to be linear or bilinear or hordijk")

        this = new(E, nu, ft, mu, wc, ws, softcurve)
        return this
    end
end

# Returns the element type that works with this material model
matching_elem_type(::SmearedCrack) = MechSolid

# Type of corresponding state structure
ip_state_type(mat::SmearedCrack) = SmearedCrackIpState


function calc_ςmax(mat::SmearedCrack, ipd::SmearedCrackIpState, upa::Float64)
    wc  = mat.wc
    ws  = mat.ws
    if mat.softcurve == "linear"
        if upa < wc
            a = mat.ft
            b = mat.ft/wc
        else
            a = b = 0.0
        end
        ςmax = a - b*upa
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
        ςmax = a - b*upa
    elseif mat.softcurve == "hordijk"
        if upa < wc
            z = (1 + 27*(upa/wc)^3)*exp(-6.93*upa/wc) - 28*(upa/wc)*exp(-6.93)
        else
            z = 0.0
        end
        ςmax = z*mat.ft
    end
    return ςmax
end


function ςmax_deriv(mat::SmearedCrack, ipd::SmearedCrackIpState, upa::Float64)
    # ∂ςmax/∂upa = dςmax
    wc  = mat.wc
    ws  = mat.ws
    if mat.softcurve == "linear"
        if upa < wc
            b = mat.ft/wc
        else
            b = 0.0
        end
        dςmax = -b
    elseif mat.softcurve == "bilinear"
        σs = 0.25*mat.ft
        if upa < ws
            b  = (mat.ft - σs)/ws
        elseif upa < wc
            b  = σs/(wc-ws)
        else
            b = 0.0
        end
        dςmax = -b
    elseif mat.softcurve == "hordijk"
        if upa < wc
            dz = ((81*upa^2*exp(-6.93*upa/wc)/wc^3) - (6.93*(1 + 27*upa^3/wc^3)*exp(-6.93*upa/wc)/wc) - 0.02738402432/wc)
        else
            dz = 0.0
        end
        dςmax = dz*mat.ft
    end
    return dςmax
end


function yield_func(mat::SmearedCrack, ipd::SmearedCrackIpState, ς::Array{Float64,1})
    ςmax = calc_ςmax(mat, ipd, ipd.upa)
    return sqrt(ς[2]^2 + ς[3]^2) + (ς[1]-ςmax)*mat.μ
end


function yield_deriv(mat::SmearedCrack, ipd::SmearedCrackIpState, ς::Array{Float64,1})
    τ = sqrt(ς[2]^2 + ς[3]^2)
    τ==0.0 && return [ mat.μ, 0.0, 0.0 ]
    return [ mat.μ, ς[2]/τ, ς[3]/τ]
end


function potential_derivs(mat::SmearedCrack, ipd::SmearedCrackIpState, ς::Array{Float64,1})
    if ς[1] >= 0.0 
        # G1:
        r = [ 2.0*ς[1]*mat.μ^2, 2.0*ς[2], 2.0*ς[3]]
    else
        # G2:
        r = [ 0.0, 2.0*ς[2], 2.0*ς[3] ]
    end
    return r
end


function calc_nu(mat::SmearedCrack, ipd::SmearedCrackIpState)
    wc  = mat.wc
    upa = ipd.upa
    if upa < wc
        z = (1 + 27*(upa/wc)^3)*exp(-6.93*upa/wc) - 28*(upa/wc)*exp(-6.93)
    else
        z = 0.0
    end
    return mat.ν*z
end

function calc_kn_ks(mat::SmearedCrack, ipd::SmearedCrackIpState)
    #kn = mat.E/ipd.h
    ν = mat.ν
    #ν = ipd.upa==0 ? mat.ν : 0.0
    ν = calc_nu(mat, ipd)
    kn = mat.E*(1-mat.ν)/((1+mat.ν)*(1-2*mat.ν))/ipd.h

    G  = mat.E/(2.0*(1.0+mat.ν))
    ks = G/ipd.h

    return kn, ks
end


function calc_Δλ(mat::SmearedCrack, ipd::SmearedCrackIpState, ςtr::Array{Float64,1})
    maxits = 100
    Δλ     = 0.0
    f      = 0.0
    upa    = 0.0
    tol    = 1e-4
    ς      = zeros(3)

    for i=1:maxits
        μ      = mat.μ
        kn, ks = calc_kn_ks(mat, ipd)

        # quantities at n+1
        if ςtr[1]>0
             ς     = [ ςtr[1]/(1+2*Δλ*kn*μ^2),  ςtr[2]/(1+2*Δλ*ks),  ςtr[3]/(1+2*Δλ*ks) ]
             dσdΔλ = [ -2*kn*μ^2*ςtr[1]/(1+2*Δλ*kn*μ^2)^2,  -2*ks*ςtr[2]/(1+2*Δλ*ks)^2,  -2*ks*ςtr[3]/(1+2*Δλ*ks)^2 ]
             drdΔλ = [ -4*kn*μ^4*ςtr[1]/(1+2*Δλ*kn*μ^2)^2,  -4*ks*ςtr[2]/(1+2*Δλ*ks)^2,  -4*ks*ςtr[3]/(1+2*Δλ*ks)^2 ]
        else
             ς     = [ ςtr[1],  ςtr[2]/(1+2*Δλ*ks),  ςtr[3]/(1+2*Δλ*ks) ]
             dσdΔλ = [ 0,  -2*ks*ςtr[2]/(1+2*Δλ*ks)^2,  -2*ks*ςtr[3]/(1+2*Δλ*ks)^2 ]
             drdΔλ = [ 0,  -4*ks*ςtr[2]/(1+2*Δλ*ks)^2,  -4*ks*ςtr[3]/(1+2*Δλ*ks)^2 ]
        end

        r      = potential_derivs(mat, ipd, ς)
        norm_r = norm(r)
        upa    = ipd.upa + Δλ*norm_r
        ςmax   = calc_ςmax(mat, ipd, upa)
        m      = ςmax_deriv(mat, ipd, upa)
        dςmaxdΔλ = m*(norm_r + Δλ*dot(r/norm_r, drdΔλ))

        f = sqrt(ς[2]^2 + ς[3]^2) + (ς[1]-ςmax)*μ
        if (ς[2]==0 && ς[3]==0) 
            dfdΔλ = (dσdΔλ[1] - dςmaxdΔλ)*μ
        else
            dfdΔλ = 1/sqrt(ς[2]^2 + ς[3]^2) * (ς[2]*dσdΔλ[2] + ς[3]*dσdΔλ[3]) + (dσdΔλ[1] - dςmaxdΔλ)*μ
        end

        Δλ = Δλ - f/dfdΔλ

        abs(f) < tol && break

        if i == maxits || isnan(Δλ)
            return Δλ, ς, upa, failure("SmearedCrack: Could not find Δλ")
        end
    end
    return Δλ, ς, upa, success()
end

function calcDe(mat::SmearedCrack, ipd::SmearedCrackIpState)
    E = mat.E
    ν = mat.ν
    #ν = ipd.upa==0 ? mat.ν : 0.0
    ν = calc_nu(mat, ipd)
    De = calcDe(E, ν, ipd.env.modeltype)

    return De
end

function calcD(mat::SmearedCrack, ipd::SmearedCrackIpState)
    h = ipd.h
    D = calcDe(mat, ipd)

    # If has no crack return complete elastic tensor
    ipd.upa==0 && return D

    ςmax = calc_ςmax(mat, ipd, ipd.upa)
    kn, ks = calc_kn_ks(mat, ipd)

    if ipd.Δλ == 0.0 # Cracked but in elastic regime
        #@show "elastic"

        D = [ kn*h  D[1,2] D[1,3] 0.0    0.0    0.0
              D[2,1] D[2,2] D[2,3] 0.0    0.0    0.0 
              D[3,1] D[3,2] D[3,3] 0.0    0.0    0.0
              0.0  0.0    0.0    D[4,4] 0.0    0.0
              0.0  0.0    0.0    0.0    2*ks*h 0.0
              0.0  0.0    0.0    0.0    0.0    2*ks*h ]
    else
        #if ςmax == 0.0 && ipd.w[1] >= 0.0
        if ςmax == 0.0
            #@show "fully cracked"
            return D*1e-5
            #z = kn*ipd.h*1e-10

            #D = [ z    0.0    0.0    0.0    0.0    0.0
                  #0.0  D[2,2] D[2,3] 0.0    0.0    0.0 
                  #0.0  D[3,2] D[3,3] 0.0    0.0    0.0
                  #0.0  0.0    0.0    D[4,4] 0.0    0.0
                  #0.0  0.0    0.0    0.0    2*z    0.0
                  #0.0  0.0    0.0    0.0    0.0    2*z ]

        else # Elastoplastic tensor D
            #@show "partially cracked"
            idx = [1, 5, 6]          # indexes for components in the crack plane
            M   = (1.0, SR2, SR2)  # Mandel's notation correction coefficients
            ς = (ipd.T*ipd.σ)[idx]./M  # stress components in the crack plane

            v = yield_deriv(mat, ipd, ς)
            r = potential_derivs(mat, ipd, ς)
            y = -mat.μ # ∂F/∂ςmax
            m = ςmax_deriv(mat, ipd, ipd.upa)  # ∂ςmax/∂upa


            den = kn*r[1]*v[1] + ks*r[2]*v[2] + ks*r[3]*v[3] - y*m*norm(r)
            Dcr = ipd.h*[  kn-kn^2*r[1]*v[1]/den    -kn*ks*r[1]*v[2]/den    -kn*ks*r[1]*v[3]/den
                          -kn*ks*r[2]*v[1]/den       ks-ks^2*r[2]*v[2]/den  -ks^2*r[2]*v[3]/den
                          -kn*ks*r[3]*v[1]/den      -ks^2*r[3]*v[2]/den      ks-ks^2*r[3]*v[3]/den ]

            D = [ Dcr[1,1]     0.0    0.0    0.0    Dcr[1,2]*SR2 Dcr[1,3]*SR2
                  #D[2,1] D[2,2] D[2,3] 0.0    0.0    0.0 
                  #D[3,1] D[3,2] D[3,3] 0.0    0.0    0.0
                  #0.0  0.0    0.0    D[4,4] 0.0    0.0
                  #0.0  kn*h   0.0    0.0    0.0    0.0 
                  #0.0  0.0    kn*h   0.0    0.0    0.0  
                  #0.0  0.0    0.0    2*ks*h  0.0    0.0
                  0.0          D[2,2] D[2,3] 0.0    0.0          0.0 
                  0.0          D[3,2] D[3,3] 0.0    0.0          0.0
                  0.0          0.0    0.0    D[4,4] 0.0          0.0
                  Dcr[2,1]*SR2 0.0    0.0    0.0    Dcr[2,2]*2   Dcr[2,3]*2
                  Dcr[3,1]*SR2 0.0    0.0    0.0    Dcr[3,2]*2   Dcr[3,3]*2 ]

        end
    end

    # Rotate D to the global xyz system
    D = ipd.T'*D*ipd.T
    return D
end


function stress_update(mat::SmearedCrack, ipd::SmearedCrackIpState, Δε::Array{Float64,1})
    # σ : stress tensor in the xyz system
    # ε : strain tensor in the xyz system
    # ς : stress vector at crack plane

    h = ipd.h
    kn, ks = calc_kn_ks(mat, ipd)
    M = (1.0, SR2, SR2)
    idx  = [1, 5, 6]     # indexes for components in the crack plane

    σini = ipd.σ
    D = calcDe(mat, ipd)
    Δσtr = D*Δε
    σtr  = ipd.σ + Δσtr
    ipd.ε += Δε

    # Check for cracks
    if ipd.upa==0
        P, V = eigen(σtr)
        if P[1]>mat.ft
            ipd.T = tensor_rot(V)
            ipd.upa = 1e-20
        else
            ipd.σ = σtr
            return Δσtr, success()
        end
    end

    T = ipd.T
    ςtr  = (T*σtr)[idx]./M
    Δw   = (T*Δε)[idx]./M*ipd.h
    ipd.w += Δw
    Ftr  = yield_func(mat, ipd, ςtr)
    ςmax = calc_ςmax(mat, ipd, ipd.upa)

    w1 = (ipd.T*ipd.ε)[1]*ipd.h

    #if ςmax == 0.0 && ipd.w[1] >= 0.0
    if ςmax == 0.0 && w1>=0
        #@show "fully cracked"
        r1 = [ ςtr[1]/kn, ςtr[2]/ks, ςtr[3]/ks ]
        ipd.Δλ   = norm(r1)
        r = normalize!(r1)
        ipd.upa += ipd.Δλ
        Decr = [kn, ks, ks]
        ς = ςtr - ipd.Δλ*Decr.*r
        if norm(ς)>1e-2
            @show ς
        end
        σtr = zeros(6)
    elseif Ftr <= 0
        #@show "elastic"
        ipd.Δλ = 0.0
        ς = ςtr
    else
        #@show "partially cracked"
        ipd.Δλ, ς, ipd.upa, status = calc_Δλ(mat, ipd, ςtr) 
        failed(status) && return ς, status
    end

    σcr = T*σtr
    σcr[idx] .= ς.*M

    wc  = mat.wc
    upa = ipd.upa
    if upa < wc
        z = (1 + 27*(upa/wc)^3)*exp(-6.93*upa/wc) - 28*(upa/wc)*exp(-6.93)
        dz = ((81*upa^2*exp(-6.93*upa/wc)/wc^3) - (6.93*(1 + 27*upa^3/wc^3)*exp(-6.93*upa/wc)/wc) - 0.02738402432/wc)
    else
        z = 0.0
        dz = 0.0
    end

    #σcr[2] *= z
    #σcr[3] *= z
    #σcr[4] *= z

    ipd.σ = T'*σcr
    Δσ    = ipd.σ - σini

    return Δσ, success()
end

function ip_state_vals(mat::SmearedCrack, ipd::SmearedCrackIpState)
    ndim  = ipd.env.ndim
    σ, ε  = ipd.σ, ipd.ε

    D = stress_strain_dict(σ, ε, ipd.env.modeltype)
    mand = (1.0, SR2, SR2)
    T    = ipd.T
    if ipd.upa>0
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
    D[:crack]  = float(ipd.upa>0)
    return D
end
