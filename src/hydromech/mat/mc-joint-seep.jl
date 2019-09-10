 #This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export MCJointSeep

mutable struct MCJointSeepIpState<:IpState
    env::ModelEnv
    σ   ::Array{Float64,1}  # stress
    w   ::Array{Float64,1}  # relative displacements
    V   ::Array{Float64,1}  # fluid velocity
    D   ::Array{Float64,1}  # distance traveled by the fluid
    L   ::Array{Float64,1} 
    S   ::Array{Float64,1}
    uw  ::Array{Float64,1}  # interface pore pressure
    upa ::Float64           # effective plastic relative displacement
    Δλ  ::Float64           # plastic multiplier
    h   ::Float64           # characteristic length from bulk elements
    function MCJointSeepIpState(env::ModelEnv=ModelEnv())
        this = new(env)
        ndim = env.ndim
        this.σ   = zeros(ndim)
        this.w   = zeros(ndim)
        this.V   = zeros(2) 
        this.D   = zeros(2) 
        this.L   = zeros(ndim-1)
        this.S   = zeros(ndim-1)
        this.uw  = zeros(3) 
        this.upa = 0.0
        this.Δλ  = 0.0
        this.h   = 0.0
        return this
    end
end

mutable struct MCJointSeep<:Material
    E  ::Float64       # Young's modulus
    nu ::Float64       # Poisson ratio
    σmax0::Float64     # tensile strength (internal variable)
    μ  ::Float64       # tangent of friction angle
    ζ  ::Float64       # factor ζ controls the elastic relative displacements 
    wc ::Float64       # critical crack opening
    ws ::Float64       # openning at inflection (where the curve slope changes)
    softcurve::String  # softening curve model ("linear" or bilinear" or "hordijk")
    flaw::Bool         # interface element break (true or false)
    k  ::Float64       # specific permeability
    γw ::Float64       # specific weight of the fluid
    α  ::Float64       # Biot's coefficient
    S  ::Float64       # Storativity coefficient
    β  ::Float64       # compressibility of fluid
    η  ::Float64       # viscosity
    kt ::Float64       # transverse leak-off coefficient
    kl ::Float64       # initial fracture opening (longitudinal flow)

    function MCJointSeep(prms::Dict{Symbol,Float64})
        return  MCJointSeep(;prms...)
    end

     function MCJointSeep(;E=NaN, nu=NaN, ft=NaN, mu=NaN, zeta=NaN, wc=NaN, ws=NaN, GF=NaN, Gf=NaN, softcurve="bilinear", flaw=false, k=NaN, kappa=NaN, gammaw=NaN, alpha=1.0, S=0.0, n=NaN, Ks=NaN, Kw=NaN, beta=0.0, eta=NaN, kt=NaN, kl=0.0)  

        !(isnan(GF) || GF>0) && error("Invalid value for GF: $GF")
        !(isnan(Gf) || Gf>0) && error("Invalid value for Gf: $Gf")

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

        !(isnan(kappa) || kappa>0) && error("Invalid value for kappa: $kappa")

        if isnan(k) 
            k = (kappa*gammaw)/eta # specific permeability = (intrinsic permeability * fluid specific weight)/viscosity
        end

        if isnan(S) 
            S = (alpha - n)/Ks + n/Kw # S = (alpha - porosity)/(bulk module of the solid) + (porosity)/(bulk module of the fluid) 
        end

        E>0.0       || error("Invalid value for E: $E")
        0<=nu<0.5   || error("Invalid value for nu: $nu")
        ft>=0       || error("Invalid value for ft: $ft")
        mu>0        || error("Invalid value for mu: $mu")
        zeta>0      || error("Invalid value for zeta: $zeta")
        wc>0        || error("Invalid value for wc: $wc")
        (isnan(ws)  || ws>0) || error("Invalid value for ws: $ws")
        (softcurve=="linear" || softcurve=="bilinear" || softcurve=="hordijk") || error("Invalid softcurve: softcurve must to be linear or bilinear or hordijk")
        k>0         || error("Invalid value for k: $k")
        gammaw>0    || error("Invalid value for gammaw: $gammaw")
        0<alpha<=1.0|| error("Invalid value for alpha: $alpha")
        S>=0.0      || error("Invalid value for S: $S")
        beta>=0     || error("Invalid value for beta: $beta")
        eta>=0      || error("Invalid value for eta: $eta")
        kt>=0       || error("Invalid value for kt: $kt")
        kl>=0       || error("Invalid value for kl: $kl")

        this = new(E, nu, ft, mu, zeta, wc, ws, softcurve, flaw, k, gammaw, alpha, S, beta, eta, kt, kl)
        return this
    end
end

# Returns the element type that works with this material model
matching_elem_type(::MCJointSeep) = HydroMechJoint

# Type of corresponding state structure
ip_state_type(mat::MCJointSeep) = MCJointSeepIpState


function yield_func(mat::MCJointSeep, ipd::MCJointSeepIpState, σ::Array{Float64,1})
    ndim = ipd.env.ndim
    σmax = calc_σmax(mat, ipd, ipd.upa)
    if ndim == 3
        return sqrt(σ[2]^2 + σ[3]^2) + (σ[1]-σmax)*mat.μ
    else
        return abs(σ[2]) + (σ[1]-σmax)*mat.μ
    end
end


function yield_deriv(mat::MCJointSeep, ipd::MCJointSeepIpState)
    ndim = ipd.env.ndim
    if ndim == 3
        return [ mat.μ, ipd.σ[2]/sqrt(ipd.σ[2]^2 + ipd.σ[3]^2), ipd.σ[3]/sqrt(ipd.σ[2]^2 + ipd.σ[3]^2)]
    else
        return [ mat.μ, sign(ipd.σ[2]) ]
    end
end


function potential_derivs(mat::MCJointSeep, ipd::MCJointSeepIpState, σ::Array{Float64,1})
    ndim = ipd.env.ndim
    if ndim == 3
            if σ[1] >= 0.0 
                # G1:
                r = [ 2.0*σ[1]*mat.μ^2, 2.0*σ[2], 2.0*σ[3]]
            else
                # G2:
                r = [ 0.0, 2.0*σ[2], 2.0*σ[3] ]
            end
    else
            if σ[1] >= 0.0 
                # G1:
                r = [ 2*σ[1]*mat.μ^2, 2*σ[2]]
            else
                # G2:
                r = [ 0.0, 2*σ[2] ]
            end
    end
    return r
end


function calc_σmax(mat::MCJointSeep, ipd::MCJointSeepIpState, upa::Float64)
	if mat.softcurve == "linear"
		if upa < mat.wc
            a = mat.σmax0
            b = mat.σmax0/mat.wc
		else
            a = 0.0
            b = 0.0
        end
        σmax = a - b*upa
    elseif mat.softcurve == "bilinear"
        σs = 0.25*mat.σmax0
        if upa < mat.ws
            a  = mat.σmax0 
            b  = (mat.σmax0 - σs)/mat.ws
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
        σmax = z*mat.σmax0
    end
    return σmax
end


function σmax_deriv(mat::MCJointSeep, ipd::MCJointSeepIpState, upa::Float64)
   # ∂σmax/∂upa = dσmax
    if mat.softcurve == "linear"
		if upa < mat.wc
            b = mat.σmax0/mat.wc
		else
            b = 0.0
        end
        dσmax = -b
    elseif mat.softcurve == "bilinear"
        σs = 0.25*mat.σmax0
        if upa < mat.ws
            b  = (mat.σmax0 - σs)/mat.ws
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
        dσmax = dz*mat.σmax0
    end
    return dσmax
end


function calc_kn_ks_De(mat::MCJointSeep, ipd::MCJointSeepIpState)
    ndim = ipd.env.ndim
    kn = mat.E*mat.ζ/ipd.h
    G  = mat.E/(2.0*(1.0+mat.nu))
    ks = G*mat.ζ/ipd.h

    if ndim == 3
        De = [  kn  0.0  0.0
               0.0   ks  0.0
               0.0  0.0   ks ]
    else
        De = [  kn   0.0
                 0.0  ks  ]
    end

    return kn, ks, De
end


function calc_Δλ(mat::MCJointSeep, ipd::MCJointSeepIpState, σtr::Array{Float64,1})
    ndim = ipd.env.ndim
    maxits = 50
    Δλ     = 0.0
    f      = 0.0
    upa    = 0.0
    tol    = 1e-4      

    for i=1:maxits
        μ      = mat.μ
    	kn, ks, De = calc_kn_ks_De(mat, ipd)

		# quantities at n+1
		if ndim == 3
			if σtr[1]>0
			     σ     = [ σtr[1]/(1+2*Δλ*kn*μ^2),  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
			     dσdΔλ = [ -2*kn*μ^2*σtr[1]/(1+2*Δλ*kn*μ^2)^2,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2,  -2*ks*σtr[3]/(1+2*Δλ*ks)^2 ]
			     drdΔλ = [ -4*kn*μ^4*σtr[1]/(1+2*Δλ*kn*μ^2)^2,  -4*ks*σtr[2]/(1+2*Δλ*ks)^2,  -4*ks*σtr[3]/(1+2*Δλ*ks)^2 ]
			else
			     σ     = [ σtr[1],  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
			     dσdΔλ = [ 0,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2,  -2*ks*σtr[3]/(1+2*Δλ*ks)^2 ]
			     drdΔλ = [ 0,  -4*ks*σtr[2]/(1+2*Δλ*ks)^2,  -4*ks*σtr[3]/(1+2*Δλ*ks)^2 ]
			end
		else
			if σtr[1]>0
			     σ     = [ σtr[1]/(1+2*Δλ*kn*μ^2),  σtr[2]/(1+2*Δλ*ks) ]
			     dσdΔλ = [ -2*kn*μ^2*σtr[1]/(1+2*Δλ*kn*μ^2)^2,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2 ]
			     drdΔλ = [ -4*kn*μ^4*σtr[1]/(1+2*Δλ*kn*μ^2)^2,  -4*ks*σtr[2]/(1+2*Δλ*ks)^2 ]
			else
			     σ     = [ σtr[1],  σtr[2]/(1+2*Δλ*ks) ]
			     dσdΔλ = [ 0,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2 ]
			     drdΔλ = [ 0,  -4*ks*σtr[2]/(1+2*Δλ*ks)^2 ]
			 end
		end
			 	
		 r      = potential_derivs(mat, ipd, σ)
		 norm_r = norm(r)
		 upa    = ipd.upa + Δλ*norm_r
		 σmax   = calc_σmax(mat, ipd, upa)
		 m      = σmax_deriv(mat, ipd, upa)
		 dσmaxdΔλ = m*(norm_r + Δλ*dot(r/norm_r, drdΔλ))

		if ndim == 3
		    f = sqrt(σ[2]^2 + σ[3]^2) + (σ[1]-σmax)*μ
		    if (σ[2]==0 && σ[3]==0) 
		        dfdΔλ = (dσdΔλ[1] - dσmaxdΔλ)*μ		      
		    else
		        dfdΔλ = 1/sqrt(σ[2]^2 + σ[3]^2) * (σ[2]*dσdΔλ[2] + σ[3]*dσdΔλ[3]) + (dσdΔλ[1] - dσmaxdΔλ)*μ
		    end
		else
			f = abs(σ[2]) + (σ[1]-σmax)*mat.μ
			dfdΔλ = sign(σ[2])*dσdΔλ[2] + (dσdΔλ[1] - dσmaxdΔλ)*μ
		end
        
        Δλ = Δλ - f/dfdΔλ

        abs(f) < tol && break

        if i == maxits || isnan(Δλ)
            @error """MCJointSeep: Could not find Δλ. This may happen when the system
            becomes hypostatic and thus the global stiffness matrix is near syngular.
            Increasing the mesh refinement may result in a nonsingular matrix.
            """ iterations=i Δλ
            #error()
        end
    end
    return Δλ
end


function calc_σ_upa(mat::MCJointSeep, ipd::MCJointSeepIpState, σtr::Array{Float64,1})
    ndim = ipd.env.ndim
    μ = mat.μ
    kn, ks, De = calc_kn_ks_De(mat, ipd)

    if ndim == 3
        if σtr[1] > 0
            σ = [σtr[1]/(1 + 2*ipd.Δλ*kn*(μ^2)), σtr[2]/(1 + 2*ipd.Δλ*ks), σtr[3]/(1 + 2*ipd.Δλ*ks)]
        else
            σ = [σtr[1], σtr[2]/(1 + 2*ipd.Δλ*ks), σtr[3]/(1 + 2*ipd.Δλ*ks)]
        end    
    else
        if σtr[1] > 0
            σ = [σtr[1]/(1 + 2*ipd.Δλ*kn*(μ^2)), σtr[2]/(1 + 2*ipd.Δλ*ks)]
        else
            σ = [σtr[1], σtr[2]/(1 + 2*ipd.Δλ*ks)]
        end    
    end
    ipd.σ = σ
    r = potential_derivs(mat, ipd, ipd.σ)
    ipd.upa += ipd.Δλ*norm(r)
    return ipd.σ, ipd.upa
end


function mountD(mat::MCJointSeep, ipd::MCJointSeepIpState)
    if mat.flaw==true && ipd.upa < mat.wc
        ipd.upa = mat.wc
    end

    ndim = ipd.env.ndim
    kn, ks, De = calc_kn_ks_De(mat, ipd)
    σmax = calc_σmax(mat, ipd, ipd.upa)

    if ipd.Δλ == 0.0  # Elastic 
        return De
    elseif σmax == 0.0 
        Dep  = De*1e-10 
        return Dep
    else
        v    = yield_deriv(mat, ipd)
        r    = potential_derivs(mat, ipd, ipd.σ)
        y    = -mat.μ # ∂F/∂σmax
        m    = σmax_deriv(mat, ipd, ipd.upa)  # ∂σmax/∂upa

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

    	return Dep
    end
end


function stress_update(mat::MCJointSeep, ipd::MCJointSeepIpState, Δw::Array{Float64,1}, Δuw::Array{Float64,1},  G::Array{Float64,1}, BfUw::Array{Float64,1}, Δt::Float64)
    ndim = ipd.env.ndim
    σini = copy(ipd.σ)

    μ = mat.μ
    kn, ks, De = calc_kn_ks_De(mat, ipd)
    σmax = calc_σmax(mat, ipd, ipd.upa) 

    if isnan(Δw[1]) || isnan(Δw[2])
        @warn "mc_joint_seep!: Invalid value for joint displacement: Δw = $Δw"
    end

    # σ trial and F trial
    σtr  = ipd.σ + De*Δw

    Ftr  = yield_func(mat, ipd, σtr) 

    # Elastic and EP integration
    if σmax == 0.0 && ipd.w[1] >= 0.0
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
        # Pure elastic increment
        ipd.Δλ = 0.0
        ipd.σ  = copy(σtr) 

    else    
        ipd.Δλ = calc_Δλ(mat, ipd, σtr) 
        ipd.σ, ipd.upa = calc_σ_upa(mat, ipd, σtr)
                      
        # Return to surface:
        F  = yield_func(mat, ipd, ipd.σ)   
        if F > 1e-3
            @warn "stress_update: The value of the yield function is $F"
        end
    end

    ipd.w += Δw
    Δσ = ipd.σ - σini

    ipd.uw += Δuw
    ipd.V  = -mat.kt*G
    ipd.D +=  ipd.V*Δt

    # compute crack aperture
    if mat.kl == 0.0
        if ipd.upa == 0.0 || ipd.w[1] <= 0.0
            kl = 0.0
        else
            kl = ipd.w[1]
        end
    else
        if mat.kl >= ipd.w[1]
            kl = mat.kl
        else 
            kl = ipd.w[1]
        end
    end 

    ipd.L  =  ((kl^3)/(12*mat.η))*BfUw
    ipd.S +=  ipd.L*Δt

    return Δσ, ipd.V, ipd.L
end


function ip_state_vals(mat::MCJointSeep, ipd::MCJointSeepIpState)
    ndim = ipd.env.ndim
    if ndim == 3
       return OrderedDict(
          :w1  => ipd.w[1] ,
          :w2  => ipd.w[2] ,
          :w3  => ipd.w[3] ,
          :s1  => ipd.σ[1] ,
          :s2  => ipd.σ[2] ,
          :s3  => ipd.σ[3] ,
          :upa => ipd.upa
          )
    else
        return OrderedDict(
          :w1  => ipd.w[1] ,
          :w2  => ipd.w[2] ,
          :s1  => ipd.σ[1] ,
          :s2  => ipd.σ[2] ,
          :upa => ipd.upa
          )
    end
end
