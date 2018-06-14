# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export SmearedCrack

mutable struct SmearedCrackIpState<:IpState
    shared_data::SharedAnalysisData
    σ::Tensor2
    ε::Tensor2
    w::Array{Float64,1} # relative displacements in the crack plane
    T::Tensor4    # rotation tensor
    upa::Float64  # max plastic displacement
    #w  ::Float64  # crack opening
    Δλ ::Float64  # plastic multiplier
    h  ::Float64  # element size fraction for a integration point
    hascrack::Bool
    function SmearedCrackIpState(shared_data::SharedAnalysisData=SharedAnalysisData()) 
        this = new(shared_data)
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
    α  ::Float64  # elastic displacement scale factor
    wc ::Float64  # critical openning
    ws ::Float64  # openning at inflection
    softcurve::String # softening model ("bilinear" or "hordijk")

    function SmearedCrack(prms::Dict{Symbol,Float64})
        return SmearedCrack(;prms...)
    end

    function SmearedCrack(;E=NaN, nu=NaN, ft=NaN, mu=NaN, alpha=NaN, wc=NaN, ws=NaN, softcurve="hordijk")
        this = new(E, nu, ft, mu, alpha, wc, ws, softcurve)
        this.α = 1.0
        this.α = 0.5
        return this 
    end
end

# Returns the element type that works with this material model
matching_elem_type(::SmearedCrack) = MechSolid

# Create a new instance of Ip data
new_ip_state(mat::SmearedCrack, shared_data::SharedAnalysisData) = SmearedCrackIpState(shared_data)

function set_state(ipd::SmearedCrackIpState; sig=zeros(0), eps=zeros(0))
    if length(sig)==6
        ipd.σ[:] = sig
    else
        if length(sig)!=0; error("SmearedCrack: Wrong size for stress array: $sig") end
    end
    if length(eps)==6
        ipd.ε[:] = eps
    else
        if length(eps)!=0; error("SmearedCrack: Wrong size for strain array: $eps") end
    end
end


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
            z = (1 + 27*(upa/wc)^3)*e^(-6.93*upa/wc) - 28*(upa/wc)*e^(-6.93)
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
    τr = sqrt(ς[2]^2 + ς[3]^2)
    τr == 0.0 && return [ mat.μ, 0.0, 0.0 ]
    return [ mat.μ, ς[2]/τr, ς[3]/τr]
end


function potential_derivs(mat::SmearedCrack, ipd::SmearedCrackIpState, ς::Array{Float64,1})
    if ipd.upa < mat.wc
        if ς[1] >= 0.0 # G1:
            r = [ 2.0*ς[1]*mat.μ^2, 2.0*ς[2], 2.0*ς[3]]
        else # G2:
            r = [ 0.0, 2.0*ς[2], 2.0*ς[3] ]
        end
    else # σmax==0.0    
        r = [ 0.0, 2.0*ς[2], 2.0*ς[3] ]
    end
    return r
end


function calcD(mat::SmearedCrack, ipd::SmearedCrackIpState)
    De = calcDe(mat.E, mat.ν, ipd.shared_data.model_type)
    ipd.hascrack || return De
    #@show ipd.hascrack

    #println("\nDe\n")
    #display(De)

    kn = mat.E*mat.α/ipd.h
    ks = kn/((1.0+mat.ν)) # for conventional eng. stress/strain ks = E/(2*(1+ν))*α/h
    #ks = E/(1.0+mat.ν) # for conventional eng. stress/strain ks = E/(2*(1+ν))
    T = ipd.T
    h = ipd.h

    if ipd.Δλ == 0.0  # Elastic 
        D = De
        D[1,:] .= ( kn*h, 0.0, 0.0, 0.0, 0.0, 0.0 )
        D[5,:] .= (    0.0, 0.0, 0.0, 0.0, ks*h, 0.0 )
        D[6,:] .= (    0.0, 0.0, 0.0, 0.0, 0.0, ks*h )

        D = T'*D*T
    else
        if ipd.upa/mat.wc>0.99 && ipd.w[1] >= 0.0 #0.95*mat.wc
            #@show "hiiiiiiiiiiiiiiiiiiiiiiiii"
            k = 1e-2
            D = De
            D[1,:] .= (   k, 0.0, 0.0, 0.0, 0.0, 0.0 )
            D[5,:] .= ( 0.0, 0.0, 0.0, 0.0,   k, 0.0 )
            D[6,:] .= ( 0.0, 0.0, 0.0, 0.0, 0.0,   k )
        #D[1,:] .= ( kn*h, 0.0, 0.0, 0.0, 0.0, 0.0 )
        #D[5,:] .= (    0.0, 0.0, 0.0, 0.0, ks*h, 0.0 )
        #D[6,:] .= (    0.0, 0.0, 0.0, 0.0, 0.0, ks*h )

            D = T'*D*T
        else

            idx  = [1, 5, 6]     # indexes for components in the crack plane
            mand = (1.0, SR2, SR2)
            ς = (T*ipd.σ)[idx]./mand  # stress components in the crack plane
            #@show T
            #@show ipd.σ
            #@show (T*ipd.σ)[idx]
            #@show ς

            v = yield_deriv(mat, ipd, ς)
            r = potential_derivs(mat, ipd, ς)
            y = -mat.μ # ∂F/∂σmax
            m = σmax_deriv(mat, ipd, ipd.upa)  # ∂σmax/∂upa
            #@show r, m

            den = kn*r[1]*v[1] + ks*r[2]*v[2] + ks*r[3]*v[3] - y*m*norm(r)
            Dcr =  ipd.h*[  kn - kn^2*r[1]*v[1]/den    -kn*ks*r[1]*v[2]/den      -kn*ks*r[1]*v[3]/den
                     -kn*ks*r[2]*v[1]/den         ks - ks^2*r[2]*v[2]/den  -ks^2*r[2]*v[3]/den
                     -kn*ks*r[3]*v[1]/den        -ks^2*r[3]*v[2]/den        ks - ks^2*r[3]*v[3]/den ]
            #@show Dcr
            #@show den
            #@show m
            #@show norm(r)
            #@show r
            #@show v
            #@show ipd.upa
            #@show ks

            D = De
            D[1,:] .= ( Dcr[1,1]    , 0.0, 0.0, 0.0, Dcr[1,2]*SR2, Dcr[1,3]*SR2 )
            D[5,:] .= ( Dcr[2,1]*SR2, 0.0, 0.0, 0.0, Dcr[2,2]    , Dcr[2,3] )
            D[6,:] .= ( Dcr[3,1]*SR2, 0.0, 0.0, 0.0, Dcr[3,2]    , Dcr[3,3] )

            D = T'*D*T

            #println("\nD\n")
            #display(D)
            #exit()
        end
    end

    #println("\nD")
    #display(D)

    return D
end


function find_intersection(mat::SmearedCrack, ipd::SmearedCrackIpState, F1::Float64, F2::Float64, σ0::Array{Float64,1}, Δσ::Array{Float64,1})
    @assert(F1*F2<0.0)

    # Regula Falsi method
    α  = 0.0
    α1 = 0.0
    α2 = 1.0
    F  = F1
    maxits = 40
    for i=1:maxits
        α  = α1 + F1/(F1-F2)*(α2-α1)
        F  = yield_func(mat, ipd, σ0 + α*Δσ)
        abs(F)<1e-7 && break

        if F<0.0
            α1 = α
            F1 = F
        else
            α2 = α
            F2 = F
        end
    end

    return α
end


function yield_n1_and_deriv(mat::SmearedCrack, ipd::SmearedCrackIpState, ςtr::Array{Float64,1}, Δλ::Float64)
    kn = mat.E*mat.α/ipd.h
    ks = kn/((1.0+mat.ν)) # for conventional eng. stress/strain ks = E/(2*(1+ν))*α/h
    μ = mat.μ

    # quantities at n+1
    if ςtr[1]>0
        ς     = [ ςtr[1]/(1+2*Δλ*kn*μ^2),  ςtr[2]/(2*ks*Δλ+1),  ςtr[3]/(2*ks*Δλ+1) ]
        dςdΔλ = [ -2*kn*μ^2*ςtr[1]/(1+2*Δλ*kn*μ^2)^2,  -2*ks*ςtr[2]/(2*ks*Δλ+1)^2,  -2*ks*ςtr[3]/(2*ks*Δλ+1)^2 ]
        drdΔλ  = [ -4*kn*μ^4*ςtr[1]/(1+2*Δλ*kn*μ^2)^2,  -4*ks*ςtr[2]/(2*ks*Δλ+1)^2,  -4*ks*ςtr[3]/(2*ks*Δλ+1)^2 ]
    else
        ς     = [ ςtr[1],  ςtr[2]/(2*ks*Δλ+1),  ςtr[3]/(2*ks*Δλ+1) ]
        dςdΔλ = [ 0,  -2*ks*ςtr[2]/(2*ks*Δλ+1)^2,  -2*ks*ςtr[3]/(2*ks*Δλ+1)^2 ]
        drdΔλ  = [ 0,  -4*ks*ςtr[2]/(2*ks*Δλ+1)^2,  -4*ks*ςtr[3]/(2*ks*Δλ+1)^2 ]
    end

    r      = potential_derivs(mat, ipd, ς)
    norm_r = norm(r)

    upa  = ipd.upa + Δλ*norm_r
    σmax = calc_σmax(mat, ipd, upa)
    fn1  = sqrt(ς[2]^2 + ς[3]^2) + (ς[1]-σmax)*μ
    m    = σmax_deriv(mat, ipd, upa)
    dςmaxdΔλ = m * (norm_r + Δλ*dot(r/norm_r, drdΔλ))


    if ς[2]==0 && ς[3]==0
        dfn1 = (dςdΔλ[1] - dςmaxdΔλ)*μ
    else
        dfn1 = 1/sqrt(ς[2]^2 + ς[3]^2) * (ς[2]*dςdΔλ[2] + ς[3]*dςdΔλ[3]) + (dςdΔλ[1] - dςmaxdΔλ)*μ
    end

    #@show fn1
    #@show dfn1
    #=
    if isnan(dfn1) || isnan(upa)
        @show upa
        @show ςtr
        @show Δλ
        @show r
        @show ς
        @show m
        @show dςmaxdΔλ
        @show fn1
        @show dfn1
        #exit()
    end
    =#
    return fn1, dfn1, ς, upa
end


function stress_update(mat::SmearedCrack, ipd::SmearedCrackIpState, Δε::Array{Float64,1})
    # σ : stress tensor in the xyz system
    # ε : strain tensor in the xyz system
    # σ̄ : stress tensor in the x'y'z' system
    # ε̄ : strain tensor in the x'y'z' system
    # ς : stress vector at crack plane
    # ξ : strain vector at crack plane

    σini  = ipd.σ
    ipd.ε = ipd.ε + Δε
    De    = calcDe(mat.E, mat.ν, ipd.shared_data.model_type)
    Δσe   = De*Δε
    σtr   = ipd.σ + Δσe
    #@show De*(ipd.ε)
    #@show σtr

    # Check for cracks

    if !ipd.hascrack
        P = eigvals(σtr)
        if P[1]>mat.ft
            P, V = eig(σtr) 
            ipd.T = tensor_rot(V) 
            ipd.hascrack = true

            #mand = (1.0, SR2, SR2)
            #T = ipd.T
            #idx  = [1, 5, 6]     # indexes for components in the crack plane
            #@show σtr
            #@show eigvals(σtr)
            #@show eig(σtr)
            #@show (T*σtr)[idx]./mand
            #@show (T*σtr)
            #s = tfull(σtr)
            #@show V'*s*V
            #exit()
        end
    end

    # Update process

    if !ipd.hascrack
        ipd.σ = σtr
        return Δσe
    end

    mand = (1.0, SR2, SR2)
    T    = ipd.T
    idx  = [1, 5, 6]     # indexes for components in the crack plane
    Δw   = (T*Δε)[idx]./mand*ipd.h

    @show Δε
    @show Δw


    kn = mat.E*mat.α/ipd.h
    ks = kn/((1.0+mat.ν)) # for conventional eng. stress/strain ks = E/(2*(1+ν))*α/h

    Decr = [kn, ks, ks]
    ςini = (T*σini)[idx]./mand # initial stress components in the crack plane
    ςtr  = ςini .+ Decr.*Δw
    Ftr  = yield_func(mat, ipd, ςtr)

    h = ipd.h
    De[1,:] .= ( kn*h, 0.0, 0.0, 0.0, 0.0, 0.0 )
    De[5,:] .= (    0.0, 0.0, 0.0, 0.0, ks*h,   0.0)
    De[6,:] .= (    0.0, 0.0, 0.0, 0.0, 0.0, ks*h )

    De = T'*De*T
    Δσe   = De*Δε



    #@show Ftr

    if Ftr <= 0.0
        # Pure elastic increment
        ipd.Δλ = 0.0
        ς = ςtr
        #ξ += Δw
    else
        #@show ipd.upa/mat.wc
        #@show ipd.w[1]/mat.wc
        #@show ipd.w[1]
        #@show ςtr

        if ipd.upa/mat.wc >= 0.99 && ipd.w[1] >= 0.0 && ςtr[1]>=0.0 
            #println("111111111111111111111111111111111111")
            r1 = [ ςtr[1]/kn, ςtr[2]/ks, ςtr[3]/ks ]
            ipd.Δλ   = norm(r1)
            ipd.upa += ipd.Δλ

            r = normalize!(r1)
            #Decr = [kn, ks, ks ]
            ς = ςtr - ipd.Δλ*Decr.*r
            #ς = [0.0, 0, 0]
            #@show ς
        else
            maxits = 50
            Δλ   = 0.0
            #Δλ   = 1e-10
            upa  = ipd.upa
            tol = 1e-2

            r = potential_derivs(mat, ipd, ςini)
            m = σmax_deriv(mat, ipd, ipd.upa)  # ∂σmax/∂upa
            Δλ0 = yield_func(mat, ipd, ςtr) / (kn*r[1]-m*norm(r))
            if isnan(Δλ0)
                @show Δλ0
            end
            if Δλ0<=0
                #println("nooooooooooooooooooooooooooooooooooooooooo")
            end
            @assert Δλ0>0
            #Δλ = Δλ0
            f = 0.0

            #for Δλ in (Δλ0, 0.0, 1e-10)
            #for Δλ in (0.0, 1e-10)
            Δλ = 0.0
            for i=1:maxits
                f, df, ς, upa = yield_n1_and_deriv(mat, ipd, ςtr, Δλ)
                Δλ = Δλ - f/df
                #@show f, df, Δλ
                abs(f) < tol && break
                if i == maxits || isnan(Δλ)
                    #warn("max iterations reached")
                    #@show f, Δλ
                    Δλ = -1.0
                    break
                end
            end

            if Δλ<0 
                σmax = calc_σmax(mat, ipd, ipd.upa)
                warn("Δλ calculation failed using NR")
                @show σmax
                @show ipd.upa/mat.wc
                @show ςini
                @show Δλ
                Δλ=Δλ0

                if ςtr[1]>0
                    ς  = [ ςtr[1]/(1+2*Δλ*kn*mat.μ^2),  ςtr[2]/(2*ks*Δλ+1),  ςtr[3]/(2*ks*Δλ+1) ]
                else
                    ς  = [ ςtr[1],  ςtr[2]/(2*ks*Δλ+1),  ςtr[3]/(2*ks*Δλ+1) ]
                end
                r    = potential_derivs(mat, ipd, ς)
                upa  = ipd.upa + Δλ*norm(r)
                #@show ς
                F = yield_func(mat, ipd, ς)
                @show F
            end
            ipd.upa = upa
            ipd.Δλ  = Δλ
            #@show ς
            #@show ipd.upa
            
        

            # Return to surface:
            Fend  = yield_func(mat, ipd, ς)
            #@show Fend
            #@show Δλ
            σmax = calc_σmax(mat, ipd, ipd.upa)
            #@show upa
            #@show σmax
            #@show ς


            #if false 
            if abs(Fend) > 1e-2
                @show Fend
            end
        end
    end


    Fend  = yield_func(mat, ipd, ς)

    # get σ̅ (in the x'y'z' system)
    ε̅ = T*ipd.ε
    σ̅ = De*ε̅

    # fix using crack stresses
    σ̅[idx] .= ς.*mand

    # rotate to xyz system
    ipd.σ = T'*σ̅

    ipd.w += Δw
    Δσ = ipd.σ - σini
    #@show Δσ

    return Δσ 
end

function ip_state_vals(mat::SmearedCrack, ipd::SmearedCrackIpState)
    ndim  = ipd.shared_data.ndim
    σ, ε  = ipd.σ, ipd.ε

    D = stress_strain_dict(σ, ε, ndim)


    mand = (1.0, SR2, SR2)
    T    = ipd.T
    if ipd.hascrack
        s1    = (T*ipd.σ)[1] # trial stress components in the crack plane
        s2    = (T*ipd.σ)[5] # trial stress components in the crack plane
        s3    = (T*ipd.σ)[6] # trial stress components in the crack plane
    else
        s1 = s2 = s3 = 0.0
    end
    D[:upa] = ipd.upa
    D[:s1]  = s1
    D[:s2]  = s2
    D[:s3]  = s3
    D[:w]  = ipd.w[1]
    D[:crack]  = float(ipd.hascrack)
    return D
end
