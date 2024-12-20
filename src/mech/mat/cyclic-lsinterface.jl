# This file ips part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export CyclicRSJoint, CyclicLSJoint

mutable struct CyclicLSJointState<:IpState
    ctx::Context
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
    function CyclicLSJointState(ctx::Context)
        this         = new(ctx)
        ndim         = ctx.ndim
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

CyclicLSJoint_params = [
    FunInfo(:CyclicLSJoint, "Consitutive model for a rod-solid interface."),
    KwArgInfo(:taumax, "Shear strength", cond=:(taumax>0)),
    KwArgInfo(:taures, "Residual shear stress", cond=:(taures>=0)),
    KwArgInfo((:speak, :s1), "Peak slip", cond=:(speak>0)),
    KwArgInfo((:sres, :s3), "Residual slip", cond=:(sres>=0)),
    KwArgInfo(:alpha, "Ascending curvature parameter", 0.4, cond=:(0.0<=alpha<=1.0)),
    KwArgInfo(:beta, "Descending curvature parameter", 1.0, cond=:(0.0<=beta<=1.0)),
    KwArgInfo(:kn, "Normal stiffness", cond=:(kn>0)),
    KwArgInfo(:ks, "Shear stiffness", 0.0, cond=:(ks>=0)),
    KwArgInfo(:p, "Perimeter", cond=:(p>0)),
    ArgCond(:(taumax>taures)),
    ArgCond(:(ks>=taumax/speak)),
]

mutable struct CyclicLSJoint<:Material
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

    function CyclicLSJoint(; kwargs...)
        args = checkargs(kwargs, CyclicLSJoint_params)
        this = new(args.taumax, args.taures, args.speak, 1.1*args.speak, args.sres, args.alpha, args.beta, args.ks, args.kn, args.p)
        return this
    end

end

const CyclicRSJoint = CyclicLSJoint


# Type of corresponding state structure
compat_state_type(::Type{CyclicLSJoint}, ::Type{MechLSJoint}, ctx::Context) = CyclicLSJointState


function tau(mat::CyclicLSJoint, ips::CyclicLSJointState, s::Float64)
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


function tau_deriv(mat::CyclicLSJoint, ips::CyclicLSJointState, s::Float64)
    s = abs(s)
    s2 = ips.speak*1.1

    if s==0.0
        s1_factor = 0.01
        s = s1_factor*ips.speak   # to avoid undefined derivative
    end

    if s<=ips.speak
        return mat.α*ips.τmax/ips.speak*(s/ips.speak)^(mat.α-1)
    elseif s<s2
        return mat.ks*1e-3
    elseif s<mat.sres
        return -mat.β*(ips.τmax-ips.τres)/(mat.sres-s2)*((s-s2)/(mat.sres-s2))^(mat.β-1)
    else
        return mat.ks*1e-3
    end
end


function calcD(mat::CyclicLSJoint, ips::CyclicLSJointState)
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
        else
            dτydsy = tau_deriv(mat, ips, s)
        end
        ks = dτydsy
    end

    if ips.ctx.ndim==2
        return [  ks  0.0
                 0.0   kn ]
    else
        return  [  ks  0.0  0.0
                  0.0   kn  0.0
                  0.0  0.0   kn ]
    end
end


function update_state!(mat::CyclicLSJoint, ips::CyclicLSJointState, Δu::Vect)
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


function stress_update2(mat::CyclicLSJoint, ips::CyclicLSJointState, Δu::Vect)
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


function ip_state_vals(mat::CyclicLSJoint, ips::CyclicLSJointState)
    return OrderedDict(
      :ur   => ips.u[1] ,
      :tau  => ips.σ[1] ,
      )
end

