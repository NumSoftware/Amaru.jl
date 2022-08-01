 #This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export PJointSeep

mutable struct PJointSeepState<:IpState
    env::ModelEnv
    σ   ::Array{Float64,1}  # stress
    w   ::Array{Float64,1}  # relative displacements
    Vt  ::Array{Float64,1}  # transverse fluid velocity
    L   ::Array{Float64,1} 
    uw  ::Array{Float64,1}  # interface pore pressure
    up ::Float64           # effective plastic relative displacement
    Δλ  ::Float64           # plastic multiplier
    h   ::Float64           # characteristic length from bulk elements
    function PJointSeepState(env::ModelEnv=ModelEnv())
        this = new(env)
        ndim = env.ndim
        this.σ   = zeros(ndim)
        this.w   = zeros(ndim)
        this.Vt  = zeros(2) 
        this.L   = zeros(ndim-1)
        this.uw  = zeros(3) 
        this.up = 0.0
        this.Δλ  = 0.0
        this.h   = 0.0
        return this
    end
end

mutable struct PJointSeep<:Material
    E  ::Float64       # Young's modulus
    ν  ::Float64       # Poisson ratio
    ft::Float64        # tensile strength (internal variable)
    fc::Float64        # compressive strength (internal variable)
    ζ  ::Float64       # factor ζ controls the elastic relative displacements (formerly α)
    wc ::Float64       # critical crack opening
    ws ::Float64       # openning at inflection (where the curve slope changes)
    softcurve::String  # softening curve model ("linear" or bilinear" or "hordijk")
    γw ::Float64       # specific weight of the fluid
    β  ::Float64       # compressibility of fluid
    η  ::Float64       # viscosity
    kt ::Float64       # transverse leak-off coefficient
    w ::Float64        # initial fracture opening (longitudinal flow)
    fracture::Bool     # pre-existing fracture (true or false)

    function PJointSeep(prms::Dict{Symbol,Float64})
        return  PJointSeep(;prms...)
    end

     function PJointSeep(;E=NaN, nu=NaN, ft=NaN, fc=NaN, zeta=NaN, wc=NaN, ws=NaN, GF=NaN, Gf=NaN, softcurve="bilinear", gammaw=NaN, beta=0.0, eta=NaN, kt=NaN, w=0.0, fracture=false)  

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

        E>0.0       || error("Invalid value for E: $E")
        0<=nu<0.5   || error("Invalid value for nu: $nu")
        ft>=0       || error("Invalid value for ft: $ft")
        fc<=0       || error("Invalid value for fc: $fc")
        zeta>0      || error("Invalid value for zeta: $zeta")
        wc>0        || error("Invalid value for wc: $wc")
        (isnan(ws)  || ws>0) || error("Invalid value for ws: $ws")
        (softcurve=="linear" || softcurve=="bilinear" || softcurve=="hordijk") || error("Invalid softcurve: softcurve must to be linear or bilinear or hordijk")
        gammaw>0    || error("Invalid value for gammaw: $gammaw")
        beta>= 0    || error("Invalid value for beta: $beta")
        eta>=0      || error("Invalid value for eta: $eta")
        kt>=0       || error("Invalid value for kt: $kt")
        w>=0       || error("Invalid value for w: $w")
        (fracture==true || fracture==false) || error("Invalid fracture: fracture must to be true or false")

        this = new(E, nu, ft, fc, zeta, wc, ws, softcurve, gammaw, beta, eta, kt, w, fracture)
        return this
    end
end

# Returns the element type that works with this material model
matching_elem_type(::PJointSeep) = HydroMechJoint

# Type of corresponding state structure
ip_state_type(mat::PJointSeep) = PJointSeepState


function yield_func(mat::PJointSeep, state::PJointSeepState, σ::Array{Float64,1}, σmax::Float64)
    ndim = state.env.ndim
    fc, ft = mat.fc, mat.ft
    a = ft - √(ft^2-fc*ft)
    β = 2*a - fc

    if ndim == 3
        return β*(σ[1] - σmax) + σ[2]^2 + σ[3]^2
    else
        return β*(σ[1] - σmax) + σ[2]^2
    end
end


function yield_deriv(mat::PJointSeep, state::PJointSeepState, σ::Array{Float64,1})
    ndim = state.env.ndim
    fc, ft = mat.fc, mat.ft
    a = ft - √(ft^2-fc*ft)
    β = 2*a - fc

    if ndim == 3
        return [ β , 2*σ[2], 2*σ[3] ]
    else
        return [ β , 2*σ[2] ]
    end
end


function potential_derivs(mat::PJointSeep, state::PJointSeepState, σ::Array{Float64,1})
    ndim = state.env.ndim
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


function calc_σmax(mat::PJointSeep, state::PJointSeepState, up::Float64)
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
            e = exp(1.0)
            z = (1 + 27*(up/mat.wc)^3)*e^(-6.93*up/mat.wc) - 28*(up/mat.wc)*e^(-6.93)
        else
            z = 0.0
        end
        σmax = z*mat.ft 
    end

    return σmax
end


function deriv_σmax_upa(mat::PJointSeep, state::PJointSeepState, up::Float64)
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
            e = exp(1.0)
            dz = ((81*up^2*e^(-6.93*up/mat.wc)/mat.wc^3) - (6.93*(1 + 27*up^3/mat.wc^3)*e^(-6.93*up/mat.wc)/mat.wc) - 0.02738402432/mat.wc)
        else
            dz = 0.0
        end
        dσmax = dz*mat.ft 
    end

    return dσmax
end


function calc_kn_ks(mat::PJointSeep, state::PJointSeepState)
    kn = mat.E*mat.ζ/state.h
    G  = mat.E/(2.0*(1.0+mat.ν))
    ks = G*mat.ζ/state.h

    return kn, ks
end


function calc_Δλ(mat::PJointSeep, state::PJointSeepState, σtr::Array{Float64,1})
    ndim = state.env.ndim
    maxits = 200
    Δλ     = 0.0
    f      = 0.0
    up    = 0.0
    tol    = 1e-2
    fc, ft = mat.fc, mat.ft
    a = ft - √(ft^2-fc*ft)
    β = 2*a - fc
    nits = 0

    for i in 1:maxits
        nits = i
        kn, ks = calc_kn_ks(mat, state)

        # quantities at n+1
        if ndim == 3
            σ     = [ σtr[1]-β*kn*Δλ,  σtr[2]/(1+2*Δλ*ks), σtr[3]/(1+2*Δλ*ks) ]
            dσdΔλ = [ -2*β*kn, -2*ks*σtr[2]/(1+2*Δλ*ks)^2, -2*ks*σtr[3]/(1+2*Δλ*ks)^2 ]
            drdΔλ = [ 0, -4*ks*σtr[2]/(1+2*Δλ*ks)^2, -4*ks*σtr[3]/(1+2*Δλ*ks)^2 ]
        else
            σ     = [ σtr[1]-β*kn*Δλ, σtr[2]/(1+2*Δλ*ks)]
            dσdΔλ = [ -2*β*kn, -2*ks*σtr[2]/(1+2*Δλ*ks)^2 ]
            drdΔλ = [ 0, -4*ks*σtr[2]/(1+2*Δλ*ks)^2 ]
        end
                 
        r      = yield_deriv(mat, state, σ)
        norm_r = norm(r)
        up    = state.up + Δλ*norm_r
        σmax   = calc_σmax(mat, state, up)
        m = deriv_σmax_upa(mat, state, up)
        dσmaxdΔλ = m*(norm_r + Δλ*dot(r/norm_r, drdΔλ))

        if ndim == 3
            f = β*(σ[1] - σmax) + σ[2]^2 + σ[3]^2
            dfdσ = [ β, 2*σ[2], 2*σ[3] ]
        else
            f = β*(σ[1] - σmax) + σ[2]^2
            dfdσ = [ β, 2*σ[2] ]
        end

        
        dfdσmax = -β
        dfdΔλ = dot(dfdσ, dσdΔλ) + dfdσmax*dσmaxdΔλ
        Δλ = Δλ - f/dfdΔλ
        # @show f, Δλ
        abs(f) < tol && break


        if i == maxits || isnan(Δλ)
            warn("""PJointSeep: Could not find Δλ. This may happen when the system
            becomes hypostatic and thus the global stiffness matrix is nearly singular.
            Increasing the mesh refinement may result in a nonsingular matrix.
            """)
            warn("iterations=$i Δλ=$Δλ")
            return 0.0, failure("?")
        end
    end
    # @show nits
    # @show f, Δλ

    return Δλ, success()
end


function calc_σ_upa(mat::PJointSeep, state::PJointSeepState, σtr::Array{Float64,1})
    ndim = state.env.ndim
    kn, ks = calc_kn_ks(mat, state)
    fc, ft = mat.fc, mat.ft
    a = ft - √(ft^2-fc*ft)
    β = 2*a - fc

    if ndim == 3
        state.σ = [ σtr[1] - β*kn*state.Δλ, σtr[2]/(1 + 2*state.Δλ*ks), σtr[3]/(1 + 2*state.Δλ*ks) ]
    else
        state.σ = [ σtr[1] - β*kn*state.Δλ, σtr[2]/(1 + 2*state.Δλ*ks) ]
    end
    r = yield_deriv(mat, state, state.σ)
    state.up += state.Δλ*norm(r)
    return state.σ, state.up
end


function mountD(mat::PJointSeep, state::PJointSeepState)
    ndim = state.env.ndim
    kn, ks = calc_kn_ks(mat, state)
    σmax = calc_σmax(mat, state, state.up)

    De = diagm([kn, ks, ks][1:ndim])

    if state.Δλ == 0.0 
        return De
    elseif σmax == 0.0 
        Dep = De*1e-4
        return Dep
    else
        fc, ft = mat.fc, mat.ft
        a = ft - √(ft^2-fc*ft)
        β = 2*a - fc

        r = yield_deriv(mat, state, state.σ)
        v = r
        y = -β  # ∂F/∂σmax
        m = deriv_σmax_upa(mat, state, state.up)  # ∂σmax/∂up

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


function stress_update(mat::PJointSeep, state::PJointSeepState, Δw::Array{Float64,1}, Δuw::Array{Float64,1},  G::Array{Float64,1}, BfUw::Array{Float64,1}, Δt::Float64)
    ndim = state.env.ndim
    σini = copy(state.σ)

    kn, ks = calc_kn_ks(mat, state)
    De = diagm([kn, ks, ks][1:ndim])

    if mat.fracture 
        state.up = mat.wc
    end 
    
    σmax = calc_σmax(mat, state, state.up) 

    if isnan(Δw[1]) || isnan(Δw[2])
        alert("PJointSeep: Invalid value for joint displacement: Δw = $Δw")
    end

    # σ trial and F trial
    σtr  = state.σ + De*Δw

    Ftr  = yield_func(mat, state, σtr, σmax) 

    # Elastic and EP integration
    if σmax == 0.0 && state.w[1] >= 0.0
        if ndim==3
            r1 = [ σtr[1]/kn, σtr[2]/ks, σtr[3]/ks ]
            r = r1/norm(r1)
            state.Δλ = norm(r1)
        else
            r1 = [ σtr[1]/kn, σtr[2]/ks ]
            r = r1/norm(r1)
            state.Δλ = norm(r1)  
        end

        state.up += state.Δλ
        state.σ = σtr - state.Δλ*De*r     

    elseif Ftr <= 0.0
        # Pure elastic increment
        state.Δλ = 0.0
        state.σ  = copy(σtr) 

    else    
        state.Δλ, status = calc_Δλ(mat, state, σtr) 
        failed(status) && return state.σ, status

        state.σ, state.up = calc_σ_upa(mat, state, σtr)
                      
        # Return to surface:
        # F  = yield_func(mat, state, state.σ)   
        # F > 1e-2 && alert("PJointSeep: Yield function value ($F) outside tolerance")
    end

    state.w += Δw
    Δσ = state.σ - σini

    state.uw += Δuw
    state.Vt  = -mat.kt*G

    # compute crack aperture
    if mat.w == 0.0
        if state.up == 0.0 || state.w[1] <= 0.0 
            w = 0.0
        else
            w = state.w[1]
        end
    else
        if mat.w >= state.w[1]
            w = mat.w
        else 
            w = state.w[1]
        end
    end 

    state.L  =  ((w^3)/(12*mat.η))*BfUw

    return Δσ, state.Vt, state.L
end


function ip_state_vals(mat::PJointSeep, state::PJointSeepState)
    ndim = state.env.ndim
    if ndim == 3
       return OrderedDict(
          :w1   => state.w[1] ,
          :w2   => state.w[2] ,
          :w3   => state.w[3] ,
          :s1   => state.σ[1] ,
          :s2   => state.σ[2] ,
          :s3   => state.σ[3] ,
          :up  => state.up  ,
          :uwf  => state.uw[3],
          :vb   => state.Vt[1],
          :vt   => state.Vt[2])
    else
        return OrderedDict(
          :w1   => state.w[1] ,
          :w2   => state.w[2] ,
          :s1   => state.σ[1] ,
          :s2   => state.σ[2] ,
          :up  => state.up  ,
          :uwf  => state.uw[3],
          :vb   => state.Vt[1],
          :vt   => state.Vt[2])
    end
end
