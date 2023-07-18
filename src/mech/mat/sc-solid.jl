# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export SmearedCrack

mutable struct SmearedCrackState<:IpState
    env::ModelEnv
    σ::Vec6
    ε::Vec6
    T::Mat6x6    # rotation tensor to the crack plane
    w::Array{Float64,1} # relative displacements in the crack plane
    up::Float64  # max plastic displacement
    Δλ ::Float64  # plastic multiplier
    h  ::Float64  # element size fraction for a integration point
    function SmearedCrackState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ = zeros(6)
        this.ε = zeros(6)
        this.T = eye(6)
        this.w   = zeros(3)
        this.up = 0.0
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


# Type of corresponding state structure
ip_state_type(mat::SmearedCrack) = SmearedCrackState


function calc_ςmax(mat::SmearedCrack, state::SmearedCrackState, up::Float64)
    wc  = mat.wc
    ws  = mat.ws
    if mat.softcurve == "linear"
        if up < wc
            a = mat.ft
            b = mat.ft/wc
        else
            a = b = 0.0
        end
        ςmax = a - b*up
    elseif mat.softcurve == "bilinear"
        σs = 0.25*mat.ft
        if up < ws
            a  = mat.ft
            b  = (mat.ft - σs)/ws
        elseif up < wc
            a  = wc*σs/(wc-ws)
            b  = σs/(wc-ws)
        else
            a = b = 0.0
        end
        ςmax = a - b*up
    elseif mat.softcurve == "hordijk"
        if up < wc
            z = (1 + 27*(up/wc)^3)*exp(-6.93*up/wc) - 28*(up/wc)*exp(-6.93)
        else
            z = 0.0
        end
        ςmax = z*mat.ft
    end
    return ςmax
end


function ςmax_deriv(mat::SmearedCrack, state::SmearedCrackState, up::Float64)
    # ∂ςmax/∂up = dςmax
    wc  = mat.wc
    ws  = mat.ws
    if mat.softcurve == "linear"
        if up < wc
            b = mat.ft/wc
        else
            b = 0.0
        end
        dςmax = -b
    elseif mat.softcurve == "bilinear"
        σs = 0.25*mat.ft
        if up < ws
            b  = (mat.ft - σs)/ws
        elseif up < wc
            b  = σs/(wc-ws)
        else
            b = 0.0
        end
        dςmax = -b
    elseif mat.softcurve == "hordijk"
        if up < wc
            dz = ((81*up^2*exp(-6.93*up/wc)/wc^3) - (6.93*(1 + 27*up^3/wc^3)*exp(-6.93*up/wc)/wc) - 0.02738402432/wc)
        else
            dz = 0.0
        end
        dςmax = dz*mat.ft
    end
    return dςmax
end


function yield_func(mat::SmearedCrack, state::SmearedCrackState, ς::Array{Float64,1})
    ςmax = calc_ςmax(mat, state, state.up)
    return sqrt(ς[2]^2 + ς[3]^2) + (ς[1]-ςmax)*mat.μ
end


function yield_deriv(mat::SmearedCrack, state::SmearedCrackState, ς::Array{Float64,1})
    τ = sqrt(ς[2]^2 + ς[3]^2)
    τ==0.0 && return [ mat.μ, 0.0, 0.0 ]
    return [ mat.μ, ς[2]/τ, ς[3]/τ]
end


function potential_derivs(mat::SmearedCrack, state::SmearedCrackState, ς::Array{Float64,1})
    if ς[1] >= 0.0 
        # G1:
        r = [ 2.0*ς[1]*mat.μ^2, 2.0*ς[2], 2.0*ς[3]]
    else
        # G2:
        r = [ 0.0, 2.0*ς[2], 2.0*ς[3] ]
    end
    return r
end


function calc_nu(mat::SmearedCrack, state::SmearedCrackState)
    wc  = mat.wc
    up = state.up
    if up < wc
        z = (1 + 27*(up/wc)^3)*exp(-6.93*up/wc) - 28*(up/wc)*exp(-6.93)
    else
        z = 0.0
    end
    return mat.ν*z
end


function calc_kn_ks(mat::SmearedCrack, state::SmearedCrackState)
    #kn = mat.E/state.h
    ν = mat.ν
    #ν = state.up==0 ? mat.ν : 0.0
    ν = calc_nu(mat, state)
    kn = mat.E*(1-mat.ν)/((1+mat.ν)*(1-2*mat.ν))/state.h

    G  = mat.E/(2.0*(1.0+mat.ν))
    ks = G/state.h

    return kn, ks
end


function calc_Δλ(mat::SmearedCrack, state::SmearedCrackState, ςtr::Array{Float64,1})
    maxits = 100
    Δλ     = 0.0
    f      = 0.0
    up    = 0.0
    tol    = 1e-4
    ς      = zeros(3)

    for i in 1:maxits
        μ      = mat.μ
        kn, ks = calc_kn_ks(mat, state)

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

        r      = potential_derivs(mat, state, ς)
        norm_r = norm(r)
        up    = state.up + Δλ*norm_r
        ςmax   = calc_ςmax(mat, state, up)
        m      = ςmax_deriv(mat, state, up)
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
            return Δλ, ς, up, failure("SmearedCrack: Could not find Δλ")
        end
    end
    return Δλ, ς, up, success()
end


function calcDe(mat::SmearedCrack, state::SmearedCrackState)
    E = mat.E
    ν = mat.ν
    #ν = state.up==0 ? mat.ν : 0.0
    ν = calc_nu(mat, state)
    De = calcDe(E, ν, state.env.ana.stressmodel)

    return De
end


function calcD(mat::SmearedCrack, state::SmearedCrackState)
    h = state.h
    D = calcDe(mat, state)

    # If has no crack return complete elastic tensor
    state.up==0 && return D

    ςmax = calc_ςmax(mat, state, state.up)
    kn, ks = calc_kn_ks(mat, state)

    if state.Δλ == 0.0 # Cracked but in elastic regime
        #@show "elastic"

        D = [ kn*h  D[1,2] D[1,3] 0.0    0.0    0.0
              D[2,1] D[2,2] D[2,3] 0.0    0.0    0.0 
              D[3,1] D[3,2] D[3,3] 0.0    0.0    0.0
              0.0  0.0    0.0    D[4,4] 0.0    0.0
              0.0  0.0    0.0    0.0    2*ks*h 0.0
              0.0  0.0    0.0    0.0    0.0    2*ks*h ]
    else
        #if ςmax == 0.0 && state.w[1] >= 0.0
        if ςmax == 0.0
            #@show "fully cracked"
            return D*1e-5
            #z = kn*state.h*1e-10

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
            ς = (state.T*state.σ)[idx]./M  # stress components in the crack plane

            v = yield_deriv(mat, state, ς)
            r = potential_derivs(mat, state, ς)
            y = -mat.μ # ∂F/∂ςmax
            m = ςmax_deriv(mat, state, state.up)  # ∂ςmax/∂up


            den = kn*r[1]*v[1] + ks*r[2]*v[2] + ks*r[3]*v[3] - y*m*norm(r)
            Dcr = state.h*[  kn-kn^2*r[1]*v[1]/den    -kn*ks*r[1]*v[2]/den    -kn*ks*r[1]*v[3]/den
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
    D = state.T'*D*state.T
    return D
end


function update_state(mat::SmearedCrack, state::SmearedCrackState, Δε::Array{Float64,1})
    # σ : stress tensor in the xyz system
    # ε : strain tensor in the xyz system
    # ς : stress vector at crack plane

    h = state.h
    kn, ks = calc_kn_ks(mat, state)
    M = (1.0, SR2, SR2)
    idx  = [1, 5, 6]     # indexes for components in the crack plane

    σini = state.σ
    D = calcDe(mat, state)
    Δσtr = D*Δε
    σtr  = state.σ + Δσtr
    state.ε += Δε

    # Check for cracks
    if state.up==0
        P, V = eigen(σtr)
        if P[1]>mat.ft
            state.T = tensor_rot(V)
            state.up = 1e-20
        else
            state.σ = σtr
            return Δσtr, success()
        end
    end

    T = state.T
    ςtr  = (T*σtr)[idx]./M
    Δw   = (T*Δε)[idx]./M*state.h
    state.w += Δw
    Ftr  = yield_func(mat, state, ςtr)
    ςmax = calc_ςmax(mat, state, state.up)

    w1 = (state.T*state.ε)[1]*state.h

    #if ςmax == 0.0 && state.w[1] >= 0.0
    if ςmax == 0.0 && w1>=0
        #@show "fully cracked"
        r1 = [ ςtr[1]/kn, ςtr[2]/ks, ςtr[3]/ks ]
        state.Δλ   = norm(r1)
        r = normalize!(r1)
        state.up += state.Δλ
        Decr = [kn, ks, ks]
        ς = ςtr - state.Δλ*Decr.*r
        if norm(ς)>1e-2
            @show ς
        end
        σtr = zeros(6)
    elseif Ftr <= 0
        #@show "elastic"
        state.Δλ = 0.0
        ς = ςtr
    else
        #@show "partially cracked"
        state.Δλ, ς, state.up, status = calc_Δλ(mat, state, ςtr) 
        failed(status) && return ς, status
    end

    σcr = T*σtr
    σcr[idx] .= ς.*M

    wc  = mat.wc
    up = state.up
    if up < wc
        z = (1 + 27*(up/wc)^3)*exp(-6.93*up/wc) - 28*(up/wc)*exp(-6.93)
        dz = ((81*up^2*exp(-6.93*up/wc)/wc^3) - (6.93*(1 + 27*up^3/wc^3)*exp(-6.93*up/wc)/wc) - 0.02738402432/wc)
    else
        z = 0.0
        dz = 0.0
    end

    #σcr[2] *= z
    #σcr[3] *= z
    #σcr[4] *= z

    state.σ = T'*σcr
    Δσ    = state.σ - σini

    return Δσ, success()
end


function ip_state_vals(mat::SmearedCrack, state::SmearedCrackState)
    ndim  = state.env.ndim
    σ, ε  = state.σ, state.ε

    D = stress_strain_dict(σ, ε, state.env.ana.stressmodel)
    mand = (1.0, SR2, SR2)
    T    = state.T
    if state.up>0
        s1 = (T*state.σ)[1]
        s2 = (T*state.σ)[5]
        s3 = (T*state.σ)[6]
    else
        s1 = s2 = s3 = 0.0
    end
    D[:up] = state.up
    D[:s1]  = s1
    D[:s2]  = s2
    D[:s3]  = s3
    D[:w]   = state.w[1]
    D[:w2]  = state.w[2]
    D[:w3]  = state.w[3]
    D[:crack]  = float(state.up>0)
    return D
end
