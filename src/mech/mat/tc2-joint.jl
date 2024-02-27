 #This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export TC2Joint

mutable struct TC2JointState<:IpState
    env::ModelEnv
    σ  ::Array{Float64,1} # stress
    w  ::Array{Float64,1} # relative displacements
    up ::Float64          # effective plastic relative displacement
    Δλ ::Float64          # plastic multiplier
    h  ::Float64          # characteristic length from bulk elements
    w_rate::Float64       # w rate wrt T
    function TC2JointState(env::ModelEnv)
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

mutable struct TC2Joint<:Material
    E     ::Float64 # Young's modulus
    ν     ::Float64 # Poisson ratio
    ft    ::Float64 # tensile strength
    fc    ::Float64 # compressive strength
    ζ     ::Float64 # factor ζ controls the elastic relative displacements (formerly α)
    wc    ::Float64 # critical crack opening
    softmodel::Symbol # softening curve model (:linear or bilinear" or :hordijk or :soft)
    softcurve::Array{Float64,2} # softening curve model (:linear or bilinear" or :hordijk or :soft)
    α::Float64        # curvature coefficient
    γ::Float64        # factor for βres
    θ::Float64        # fator for the surface reduction speed
    β0::Float64       # initial shear stress for σn=0

    function TC2Joint(; args...)
        args = checkargs(args, func_params(TC2Joint))

        wc = args.wc
        softmodel = args.softmodel
        softcurve = args.softcurve

        if length(softcurve)>0
            if softcurve[end,2]==0
                wc = softcurve[end,1]
            else
                wc = Inf
            end
            softmodel = :custom
        end
        
        if wc==0.0
            GF = args.GF
            ft = args.ft
            if softmodel == :linear
                wc = round(2*GF/ft, sigdigits=5)
            elseif softmodel == :bilinear
                wc = round(5*GF/ft, sigdigits=5)
            elseif softmodel==:hordijk
                wc = round(GF/(0.1947019536*ft), sigdigits=5)  
            elseif softmodel==:soft
                wc = round(GF/(0.1947019536*ft), sigdigits=5)
            end
        end

        α = args.alpha
        ft = args.ft
        fc = args.fc

        ta = 0.05*pi
        tb = 0.5*pi

        f(t) = begin
            a = fc/2*(1-cos(t)) # negative value
            b = -fc/2*(sin(t))  # positive value
            χ = (ft-a)/ft
            β0 = b/asinh(α*χ)
            fc/2 - a + α*β0*b/(ft*√(α^2*χ^2 + 1))
        end
        
        t, _ = findroot(f, ta, tb, tol=1e-4, method=:default)
        t>pi/2 && throw(AmaruException("Invalid value for β0 was found. Check fc and ft values"))
        
        a = fc/2*(1-cos(t)) # negative value
        b = -fc/2*(sin(t))  # positive value
        χ = (ft-a)/ft
        β0 = b/asinh(α*χ)

        this = new(args.E, args.nu, args.ft, args.fc, args.zeta, wc, softmodel, softcurve, args.alpha, args.gamma, args.theta, β0)
        return this
    end
end



func_params(::Type{TC2Joint}) = [
    FunInfo( :TC2Joint, "Creates a `TC2Joint` material model.", ""),
    ArgInfo( :E, "Young modulus", condition=:(E>0)),
    ArgInfo( :nu, "Poisson ratio", condition=:(0<=nu<0.5)),
    ArgInfo( :fc, "Compressive strength", condition=:(fc<0)),
    ArgInfo( :ft, "Tensile strength", condition=:(ft>0)),
    ArgInfo( :zeta, "Joint elastic stiffness factgor", condition=:(zeta>0)),
    ArgInfo( :alpha, "Failure surface shape", 1.5, condition=:(alpha>0)),
    ArgInfo( :gamma, "Failure surface minimum size", 0.1, condition=:(gamma>=0)),
    ArgInfo( :theta, "Failure surface reduction speed", 1.5, condition=:(theta>=0)),
    ArgInfo( :wc, "Critical crack opening", 0.0, condition=:(wc>=0)),
    ArgInfo( :GF, "Fracture energy", 0.0, condition=:(GF>=0)),
    ArgInfo( :softmodel, "Softening model", :hordijk, values=(:linear, :bilinear, :hordijk, :soft, :custom), type=Symbol),
    ArgInfo( :softcurve, "Softening curve", zeros(0,0), type=Array),
]
@doc make_doc(func_params(TC2Joint)) TC2Joint


function paramsdict(mat::TC2Joint)
    mat = OrderedDict( string(field) => getfield(mat, field) for field in fieldnames(typeof(mat)) )

    if mat.softmodel in (:hordijk, :soft)
        mat["GF"] = 0.1943*mat.ft*mat.wc
    elseif mat.softmodel == :bilinear
        mat["GF"] = mat.ft*mat.wc/5
    elseif mat.softmodel == :linear
        mat["GF"] = mat.ft*mat.wc/2
    else
        mat["GF"] = NaN
    end
    return mat
end


# Type of corresponding state structure
compat_state_type(::Type{TC2Joint}, ::Type{MechJoint}, env::ModelEnv) = TC2JointState


function yield_func(mat::TC2Joint, state::TC2JointState, σ::Array{Float64,1}, σmax::Float64)
    β = calc_β(mat, σmax)
    χ = (σmax-σ[1])/mat.ft

    if state.env.ndim == 3
        τnorm = sqrt(σ[2]^2 + σ[3]^2)
    else
        τnorm = abs(σ[2])
    end

    # @show σmax
    # @show β
    # @show χ

    # @show τnorm - β*χ*exp(-mat.α*χ)
    # error()
    return τnorm - β*asinh(mat.α*χ)
end


function yield_derivs(mat::TC2Joint, state::TC2JointState, σ::Array{Float64,1}, σmax::Float64)
    ft   = mat.ft
    α    = mat.α
    β    = calc_β(mat, σmax)
    βres = mat.γ*mat.β0
    χ    = (σmax-σ[1])/mat.ft
    
    dfdσn  = α*β/(ft*√(α^2*χ^2 + 1))
    
    if state.env.ndim == 3
        τnorm = sqrt(σ[2]^2 + σ[3]^2)
        dfdσ = [ dfdσn, σ[2]/τnorm, σ[3]/τnorm]
    else
        dfdσ = [ dfdσn, sign(σ[2]) ]
    end

    if σmax>0
        dβdσmax = (β-βres)*mat.θ/σmax
    else
        dβdσmax = 0.0
    end
    dfdσmax = -dβdσmax*χ*exp(-α*χ) - β*exp(-α*χ)/ft*(1-α*χ)

    # @show dβdσmax
    # @show -dβdσmax*χ*exp(-α*χ)
    # @show β*exp(-α*χ)/ft*(1-α*χ)
    # @show dfdσmax

    return dfdσ, dfdσmax
end


function potential_derivs(mat::TC2Joint, state::TC2JointState, σ::Array{Float64,1})
    ndim = state.env.ndim
    if ndim == 3
        if σ[1] > 0.0 
            # G1:
            r = Float64[ 2.0*σ[1], 2.0*σ[2], 2.0*σ[3]]
        else
            # G2:
            r = Float64[ 0.0, 2.0*σ[2], 2.0*σ[3] ]
        end
        if r[1]==r[2]==r[3]==0.0
            r = Float64[ 1.0, 0.0, 0.0 ] # important
        end
    else
        if σ[1] > 0.0 
            # G1:
            r = Float64[ 2*σ[1], 2*σ[2]]
        else
            # G2:
            r = Float64[ 0.0, 2*σ[2] ]
        end
        if r[1]==r[2]==0.0
            r = Float64[ 1.0, 0.0 ] # important
        end
    end
    return r
end


function calc_σmax(mat::TC2Joint, state::TC2JointState, up::Float64)
    if mat.softmodel == :linear
        if up < mat.wc
            a = mat.ft 
            b = mat.ft /mat.wc
        else
            a = 0.0
            b = 0.0
        end
        σmax = a - b*up
    elseif mat.softmodel == :bilinear
        σs = 0.25*mat.ft
        ws = mat.wc*0.15
        if up < ws
            a  = mat.ft  
            b  = (mat.ft  - σs)/ws
        elseif up < mat.wc
            a  = mat.wc*σs/(mat.wc-ws)
            b  = σs/(mat.wc-ws)
        else
            a = 0.0
            b = 0.0
        end
        σmax = a - b*up
    elseif mat.softmodel == :hordijk
        if up < mat.wc
            z = (1 + 27*(up/mat.wc)^3)*exp(-6.93*up/mat.wc) - 28*(up/mat.wc)*exp(-6.93)
        else
            z = 0.0
        end
        σmax = z*mat.ft
    elseif mat.softmodel == :soft
        dσmaxdup = 0.55
        a = 1.30837
        if up == 0.0
            z = 1.0
        elseif 0.0 < up < mat.wc
            x = up/mat.wc
            z = 1.0 - a^(1.0 - 1.0/x^dσmaxdup)
        else
            z = 0.0
        end
        σmax = z*mat.ft
    elseif mat.softmodel == :custom
        σmax = interpolate(view(mat.softcurve, :, 1), view(mat.softcurve, :, 2), up)
    end

    return σmax
end


function calc_β(mat::TC2Joint, σmax::Float64)
    βres  = mat.γ*mat.β0
    return βres + (mat.β0-βres)*(σmax/mat.ft)^mat.θ
end


function deriv_σmax_up(mat::TC2Joint, state::TC2JointState, up::Float64)
    if mat.softmodel == :linear
        if up < mat.wc
            b = mat.ft /mat.wc
        else
            b = 0.0
        end
        dσmax = -b
    elseif mat.softmodel == :bilinear
        ws = mat.wc*0.15
        σs = 0.25*mat.ft 
        if up < ws
            b  = (mat.ft  - σs)/ws
        elseif up < mat.wc
            b  = σs/(mat.wc-ws)
        else
            b = 0.0
        end
        dσmax = -b
    elseif mat.softmodel == :hordijk
        if up < mat.wc
            dz = ((81*up^2*exp(-6.93*up/mat.wc)/mat.wc^3) - (6.93*(1 + 27*up^3/mat.wc^3)*exp(-6.93*up/mat.wc)/mat.wc) - 0.02738402432/mat.wc)
        else
            dz = 0.0
        end
        dσmax = dz*mat.ft 
    elseif mat.softmodel == :soft
        dσmaxdup = 0.55
        a = 1.30837

        if up == 0.0
            dz = 0.0
        elseif up < mat.wc
            x = up/mat.wc
            dz =  -dσmaxdup*log(a)*a^(1-x^-dσmaxdup)*x^(-dσmaxdup-1)/mat.wc

        else
            dz = 0.0
        end
        dσmax = dz*mat.ft 
    elseif mat.softmodel == :custom
        dσmax = derivative(view(mat.softcurve, :, 1), view(mat.softcurve, :, 2), up)
    end

    return dσmax
end


function calc_kn_ks(mat::TC2Joint, state::TC2JointState)
    kn = mat.E*mat.ζ/state.h
    G  = mat.E/(2.0*(1.0+mat.ν))
    ks = G*mat.ζ/state.h
    return kn, ks
end


function calcD(mat::TC2Joint, state::TC2JointState)

    ndim = state.env.ndim
    kn, ks = calc_kn_ks(mat, state)
    σmax = calc_σmax(mat, state, state.up)

    De = diagm([kn, ks, ks][1:ndim])

    if state.Δλ == 0.0  # Elastic 
        return De
    elseif σmax == 0.0 && state.w[1] >= 0.0
        Dep = De*1e-4
        # Dep = De*1e-3
        return Dep
    else
        r = potential_derivs(mat, state, state.σ) # ∂g/∂σ
        dfdσ, dfdσmax = yield_derivs(mat, state, state.σ, σmax)
        dσmaxdup = deriv_σmax_up(mat, state, state.up)  # ∂σmax/∂up

        if ndim == 3
            den = kn*r[1]*dfdσ[1] + ks*r[2]*dfdσ[2] + ks*r[3]*dfdσ[3] - dfdσmax*dσmaxdup*norm(r)

            Dep = [   kn - kn^2*r[1]*dfdσ[1]/den    -kn*ks*r[1]*dfdσ[2]/den      -kn*ks*r[1]*dfdσ[3]/den
                     -kn*ks*r[2]*dfdσ[1]/den         ks - ks^2*r[2]*dfdσ[2]/den  -ks^2*r[2]*dfdσ[3]/den
                     -kn*ks*r[3]*dfdσ[1]/den        -ks^2*r[3]*dfdσ[2]/den        ks - ks^2*r[3]*dfdσ[3]/den ]
        else
            den = kn*r[1]*dfdσ[1] + ks*r[2]*dfdσ[2] - dfdσmax*dσmaxdup*norm(r)

            Dep = [   kn - kn^2*r[1]*dfdσ[1]/den    -kn*ks*r[1]*dfdσ[2]/den      
                     -kn*ks*r[2]*dfdσ[1]/den         ks - ks^2*r[2]*dfdσ[2]/den  ]
        end

        return Dep
    end
end


function calc_σ_up_Δλ(mat::TC2Joint, state::TC2JointState, σtr::Array{Float64,1})
    ndim = state.env.ndim
    Δλ   = 0.0
    up   = 0.0
    σ    = zeros(ndim)
    σ0   = zeros(ndim)

    tol    = 1e-6
    maxits = 50
    for i in 1:maxits
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
        up     = state.up + Δλ*norm_r
        σmax   = calc_σmax(mat, state, up)
        f      = yield_func(mat, state, σ, σmax)
        # @show σmax
        dfdσ, dfdσmax = yield_derivs(mat, state, σ, σmax)
        # @show r
        # @show dfdσ
        # @show dfdσmax

        dσmaxdup = deriv_σmax_up(mat, state, up)
        dσmaxdΔλ = dσmaxdup*(norm_r + Δλ*dot(r/norm_r, drdΔλ))
        dfdΔλ = dot(dfdσ, dσdΔλ) + dfdσmax*dσmaxdΔλ
        Δλ = Δλ - f/dfdΔλ
        
        if Δλ<=0 || isnan(Δλ) || i==maxits
            # return 0.0, state.σ, 0.0, failure("TC2Joint: failed to find Δλ")
            # switch to bissection method
            return calc_σ_up_Δλ_bis(mat, state, σtr)
        end

        if maximum(abs, σ-σ0) <= tol
            break
        end
        σ0 .= σ
    end

    return σ, up, Δλ, success()
end


function calc_σ_up(mat::TC2Joint, state::TC2JointState, σtr::Array{Float64,1}, Δλ::Float64)
    ndim = state.env.ndim
    kn, ks  = calc_kn_ks(mat, state)

    if ndim == 3
        if σtr[1]>0
            σ = [ σtr[1]/(1+2*Δλ*kn),  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
        else
            σ = [ σtr[1],  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
        end
    else
        if σtr[1]>0
            σ = [ σtr[1]/(1+2*Δλ*kn),  σtr[2]/(1+2*Δλ*ks) ]
        else
            σ = [ σtr[1],  σtr[2]/(1+2*Δλ*ks) ]
        end
    end

    r  = potential_derivs(mat, state, σ)
    up = state.up + Δλ*norm(r)
    return σ, up
end

function calc_σ_up_Δλ_bis(mat::TC2Joint, state::TC2JointState, σtr::Array{Float64,1})
    ndim    = state.env.ndim
    kn, ks  = calc_kn_ks(mat, state)
    De      = diagm([kn, ks, ks][1:ndim])
    r       = potential_derivs(mat, state, state.σ)

    ff(Δλ)  = begin
        # quantities at n+1
        σ, up = calc_σ_up(mat, state, σtr, Δλ)
        σmax = calc_σmax(mat, state, up)
        yield_func(mat, state, σ, σmax)
    end

    # find root interval from Δλ estimative
    Δλ0 = norm(σtr-state.σ)/norm(De*r)
    a, b, status = findrootinterval(ff, 0.0, Δλ0)
    failed(status) && return state.σ, 0.0, 0.0, status

    Δλ, status = findroot(ff, a, b, ftol=1e-5, method=:bisection)
    failed(status) && return state.σ, 0.0, 0.0, status

    σ, up = calc_σ_up(mat, state, σtr, Δλ)
    return σ, up, Δλ, success()  

end


function update_state!(mat::TC2Joint, state::TC2JointState, Δw::Array{Float64,1})

    ndim = state.env.ndim
    σini = copy(state.σ)

    kn, ks = calc_kn_ks(mat, state)
    De = diagm([kn, ks, ks][1:ndim])
    σmax = calc_σmax(mat, state, state.up)  

    if isnan(Δw[1]) || isnan(Δw[2])
        alert("TC2Joint: Invalid value for joint displacement: Δw = $Δw")
    end

    # σ trial and F trial
    σtr  = state.σ + De*Δw

    Ftr  = yield_func(mat, state, σtr, σmax)

    # Elastic and EP integration
    if Ftr <= 0.0
        state.Δλ  = 0.0
        state.σ  .= σtr
    elseif state.up>=mat.wc && σtr[1]>0
        if ndim==3
            Δup = norm([ σtr[1]/kn, σtr[2]/ks, σtr[3]/ks ])
        else
            Δup = norm([ σtr[1]/kn, σtr[2]/ks ])
        end
        state.up += Δup
        state.σ  .= 0.0
        state.Δλ  = 1.0
    else
        # Plastic increment
        σ, up, Δλ, status = calc_σ_up_Δλ(mat, state, σtr)
        failed(status) && return state.σ, status

        state.σ, state.up, state.Δλ = σ, up, Δλ
    end
    state.w += Δw
    Δσ = state.σ - σini
    return Δσ, success()
end


function ip_state_vals(mat::TC2Joint, state::TC2JointState)
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


function output_keys(mat::TC2Joint)
    return Symbol[:jw1, :js1, :jup]
end