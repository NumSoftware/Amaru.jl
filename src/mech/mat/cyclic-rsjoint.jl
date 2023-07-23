# This file ips part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export CyclicRSJoint

mutable struct CyclicRSJointState<:IpState
    env    ::ModelEnv
    σ      ::Array{Float64,1}
    u      ::Array{Float64,1}
    τnl     ::Float64    # max stress
    τmax   ::Float64
    τres   ::Float64
    speak  ::Float64
    srev     ::Float64    # accumulated relative displacement
    sacum  ::Float64    # accumulated relative displacement
    sneg   ::Float64
    spos   ::Float64
    elastic::Bool
    function CyclicRSJointState()
        this         = new(env)
        ndim         = env.ndim
        this.σ       = zeros(ndim)
        this.u       = zeros(ndim)
        this.τnl      = 0.0
        this.τmax    = 0.0
        this.τres    = 0.0
        this.speak   = 0.0
        this.srev      = 0.0
        this.sacum   = 0.0
        this.sneg    = 0.0
        this.spos    = 0.0
        this.elastic = false
        return this
    end
end

mutable struct CyclicRSJoint<:Material
    τmax:: Float64
    τres:: Float64
    speak:: Float64
    s2  :: Float64
    sres:: Float64
    α   :: Float64
    β   :: Float64
    ks  :: Float64
    kn  :: Float64
    p   :: Float64

    function CyclicRSJoint(prms::Dict{Symbol,Float64})
        return  CyclicRSJoint(;prms...)
    end

    function CyclicRSJoint(;taumax=NaN, taures=NaN, speak=NaN, sres=NaN, alpha=NaN, beta=NaN, kn=NaN, ks=NaN, p=NaN, A=NaN, dm=NaN)
        
        @check speak>0

        # Estimate the perimeter p
        @check p>0 || A>0 || dm>0
        if isnan(p)
            if A>0
                p = 2.0*√(A*pi)
            else
                p = pi*dm
            end
        end

        # Estimate taumax if not provided
        @check ks>0
        isnan(taumax) && (taumax = ks*speak)
        @check taumax>=taures
        @check ks>=taumax/speak

        # Define alpha if not provided
        isnan(alpha) && (alpha = 1.0)
        @check 0.0<=alpha<=1.0

        # Define beta if not provided
        isnan(beta) && (beta = 1.0)
        @assert beta>=0.0 beta=1

        # Estimate kn if not provided
        isnan(kn) && (kn = ks)
        @check kn>0

        this = new(taumax, taures, speak, 1.1*speak, sres, alpha, beta, ks, kn, p)
        return this
    end
end

# Element types that work with this material
compat_elem_types(::Type{CyclicRSJoint}) = (MechRSJoint,)

# Type of corresponding state structure
compat_state_type(::Type{CyclicRSJoint}) = CyclicRSJointState


function tau(mat::CyclicRSJoint, ips::CyclicRSJointState, s::Float64)
    s = abs(s)
    s2 = ips.speak*1.1
    if s<ips.speak
        return ips.τmax*(s/ips.speak)^mat.α
    elseif s<s2
        return ips.τmax
    elseif s<mat.sres
        return ips.τmax - (ips.τmax-ips.τres)*((s-s2)/(mat.sres-s2))^mat.β
    else
        return ips.τres
    end
end


function tau_deriv(mat::CyclicRSJoint, ips::CyclicRSJointState, s::Float64)
    s = abs(s)
    s2 = ips.speak*1.1

    if s==0.0
        s1_factor = 0.01
        s = s1_factor*ips.speak   # to avoid undefined derivative
    end

    if s<=ips.speak
        return ips.τmax/ips.speak*(s/ips.speak)^(mat.α-1)
    elseif s<s2
        return mat.ks*1e-3
        # return 1.0
    elseif s<mat.sres
        return -(ips.τmax-ips.τres)/(mat.sres-s2)*((s-s2)/(mat.sres-s2))^(mat.β-1)
    else
        # return 1.0
        return mat.ks*1e-3
    end
end


function calcD(mat::CyclicRSJoint, ips::CyclicRSJointState)
    ndim = ips.env.ndim
    ks = mat.ks
    kn = mat.kn

    if !ips.elastic
        s = ips.u[1]
        τ = ips.σ[1]

        if ips.τmax==0.0 && ips.τres==0.0 && ips.speak==0.0
            ips.τmax  = mat.τmax
            ips.τres  = mat.τres
            ips.speak = mat.speak
        end

        if s*τ<0.0 || abs(τ)>ips.τmax
            dτydsy = 1.0    
            # @show "ks=1.0"
            # @show τ
        else
            dτydsy = tau_deriv(mat, ips, s)
        end
        ks = dτydsy
    end

    if ips.env.ndim==2
        return [  ks  0.0
                 0.0   kn ]
    else
        return  [  ks  0.0  0.0
                  0.0   kn  0.0
                  0.0  0.0   kn ]
    end
end


function update_state!(mat::CyclicRSJoint, ips::CyclicRSJointState, Δu::Vect)
    ks = mat.ks
    kn = mat.kn
    s  = ips.u[1]   # relative displacement
    Δs = Δu[1]      # relative displacement increment
    τini = ips.σ[1] # initial shear stress
    τtr  = τini + ks*Δs # elastic trial
    str  = s + Δs
    ips.sacum += abs(Δs)
    s>ips.spos && (ips.spos=s)
    s<ips.sneg && (ips.sneg=s)

    # amplitude
    sa = ips.spos-ips.sneg
    @assert sa>=0.0
    sh = sa*ips.srev/mat.speak^2

    ips.τmax = mat.τmax*exp(-0.0152*sh)

    ips.speak = mat.speak*(1 + 0.33*sh^0.5)
    
    if ips.sacum<mat.speak
        ips.τres = mat.τres*(ips.sacum/mat.speak)^mat.α
    else
        ips.τres = mat.τres*exp(-0.0027*sh)
    end


    if str*τtr<0
        τnl = ips.τres
    else
        τnl = max(ips.τres, tau(mat, ips, str))
    end
    
    ftr = abs(τtr) - τnl

    if ftr<0.0
        τ = τtr
        ips.elastic = true
    else
        if str*τtr<0
            Δsrev = (abs(τtr)-τnl)/ks
            @assert Δsrev>0.0
            ips.srev += Δsrev
        end
        τ = sign(τtr)*τnl
        Δτ = τ - τini
        ips.elastic = false
    end

    # calculate Δσ
    Δτ = τ - τini
    Δσ = kn*Δu
    Δσ[1] = Δτ

    # update u and σ
    ips.u .+= Δu
    ips.σ .+= Δσ

    return Δσ, success()
end


function stress_update2(mat::CyclicRSJoint, ips::CyclicRSJointState, Δu::Vect)
    ks = mat.ks
    kn = mat.kn
    s  = ips.u[1]   # relative displacement
    Δs = Δu[1]      # relative displacement increment
    τini = ips.σ[1] # initial shear stress
    τtr  = τini + ks*Δs # elastic trial
    str  = s + Δs
    ips.sacum += abs(Δs)
    s>ips.spos && (ips.spos=s)
    s<ips.sneg && (ips.sneg=s)

    # amplitude
    sa = ips.spos-ips.sneg
    @assert sa>=0.0
    sh = sa*ips.srev/mat.sres^2

    # τmax = mat.τmax*(1 - (min(ips.srev, mat.sres)/mat.sres)^0.8)
    # τmax = mat.τmax*exp(-1.2*(ips.srev/mat.sres))
    # τmax = mat.τmax*min(1, 1.2*exp(-1.8*(ips.spos-ips.sneg)/mat.sres))
    # ips.τmax = mat.τmax*exp(-1.25*sh)
    ips.τmax = mat.τmax*exp(-1.02*sh)

    # ips.speak = mat.speak*(1 + log(1 + 5*ips.srev/mat.sres))
    ips.speak = mat.speak*(1 + 2.8*sh^0.5)
    
    if ips.sacum<mat.speak
        ips.τres = mat.τres*(ips.sacum/mat.speak)^mat.α
    else
        ips.τres = mat.τres*exp(-0.17*sh)
    end


    if str*τtr<0
        τnl = ips.τres
    else
        τnl = max(ips.τres, tau(mat, ips, str))
    end
    
    ftr = abs(τtr) - τnl

    if ftr<0.0
        τ = τtr
        ips.elastic = true
    else
        if str*τtr<0
            Δsrev = (abs(τtr)-τnl)/ks
            @assert Δsrev>0.0
            ips.srev += Δsrev
        end
        τ = sign(τtr)*τnl
        Δτ = τ - τini
        ips.elastic = false
    end

    # calculate Δσ
    Δτ = τ - τini
    Δσ = kn*Δu
    Δσ[1] = Δτ

    # update u and σ
    ips.u .+= Δu
    ips.σ .+= Δσ

    return Δσ, success()
end


function ip_state_vals(mat::CyclicRSJoint, ips::CyclicRSJointState)
    return OrderedDict(
      :ur   => ips.u[1] ,
      :tau  => ips.σ[1] ,
      )
end

