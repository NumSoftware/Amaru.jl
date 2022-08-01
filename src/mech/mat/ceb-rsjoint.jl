# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export CebRSJoint

mutable struct CebRSJointState<:IpState
    env::ModelEnv
    σ  ::Array{Float64,1}
    u  ::Array{Float64,1}
    τy ::Float64      # max stress
    sy ::Float64      # accumulated relative displacement
    elastic::Bool
    function CebRSJointState(env::ModelEnv=ModelEnv())
        this = new(env)
        ndim = env.ndim
        this.σ = zeros(ndim)
        this.u = zeros(ndim)
        this.τy = 0.0
        this.sy = 0.0
        this.elastic = false
        return this
    end
end

mutable struct CebRSJoint<:Material
    τmax:: Float64
    τres:: Float64
    s1  :: Float64
    s2  :: Float64
    s3  :: Float64
    α   :: Float64
    β   :: Float64
    ks  :: Float64
    kn  :: Float64
    p   :: Float64

    function CebRSJoint(prms::Dict{Symbol,Float64})
        return  CebRSJoint(;prms...)
    end

    function CebRSJoint(;taumax=NaN, taures=NaN, TauM=NaN, TauR=NaN, s1=NaN, s2=NaN, s3=NaN, alpha=NaN, beta=NaN, kn=NaN, ks=NaN, p=NaN, A=NaN, dm=NaN)
        
        @check s1>0
        @check s2>s1
        @check s3>s2

        # Estimate the perimeter p
        @check p>0 || A>0 || dm>0
        if isnan(p)
            if A>0
                p = 2.0*√(A*pi)
            else
                p = pi*dm
            end
        end

        # Estimate TauM if not provided
        @check ks>0
        !isnan(TauM) && (taumax=TauM)
        !isnan(TauR) && (taures=TauR)

        isnan(taumax) && (taumax = ks*s1)
        @check taumax>=taures
        @check ks>=taumax/s1

        # Define alpha if not provided
        isnan(alpha) && (alpha = 1.0)
        @check 0.0<=alpha<=1.0

        # Define beta if not provided
        isnan(beta) && (beta = 1.0)
        @assert beta>=0.0 beta=1

        # Estimate kn if not provided
        isnan(kn) && (kn = ks)
        @check kn>0

        this = new(taumax, taures, s1, s2, s3, alpha, beta, ks, kn, p)
        return this
    end
end

# Returns the element type that works with this material
matching_elem_type(::CebRSJoint) = MechRodSolidJoint

# Creates a new instance of Ip data
ip_state_type(mat::CebRSJoint) = CebRSJointState

CEBJoint1D = CebRSJoint #! deprecated
export CEBJoint1D

function Tau(mat::CebRSJoint, sy::Float64)
    if sy<mat.s1
        return mat.τmax*(sy/mat.s1)^mat.α
    elseif sy<mat.s2
        return mat.τmax
    elseif sy<mat.s3
        return mat.τmax - (mat.τmax-mat.τres)*((sy-mat.s2)/(mat.s3-mat.s2))^mat.β
    else
        return mat.τres  + 1*(sy-mat.s3)
    end
end


function deriv(mat::CebRSJoint, ipd::CebRSJointState, sy::Float64)
    if sy==0.0
        s1_factor = 0.01
        sy = s1_factor*mat.s1   # to avoid undefined derivative
    end

    if sy<=mat.s1
        return mat.τmax/mat.s1*(sy/mat.s1)^(mat.α-1)
    elseif sy<mat.s2
        #return mat.ks/10000
        return 1.0
    elseif sy<mat.s3
        return -(mat.τmax-mat.τres)/(mat.s3-mat.s2)*((sy-mat.s2)/(mat.s3-mat.s2))^(mat.β-1)
    else
        return 1.0
        # return mat.ks/10000
    end
end


function calcD(mat::CebRSJoint, ipd::CebRSJointState)
    ndim = ipd.env.ndim
    ks = mat.ks

    if !ipd.elastic
        dτydsy = deriv(mat, ipd, ipd.sy)
        # ks = ks*dτydsy/(ks+dτydsy)
        ks = dτydsy
        # @s ks
    end

    kn = mat.kn
    if ipd.env.ndim==2
        return [  ks  0.0
                 0.0   kn ]
    else
        return  [  ks  0.0  0.0
                  0.0   kn  0.0
                  0.0  0.0   kn ]
    end
end


function yield_func(mat::CebRSJoint, ipd::CebRSJointState, τ::Float64)
    return abs(τ) - ipd.τy
end


function stress_update_n(mat::CebRSJoint, ipd::CebRSJointState, Δu::Vect)
    ks = mat.ks
    kn = mat.kn
    Δs = Δu[1]      # relative displacement
    τini = ipd.σ[1] # initial shear stress
    τtr  = τini + ks*Δs # elastic trial

    ftr  = yield_func(mat, ipd, τtr)

    if ftr<0.0
        τ = τtr
        ipd.elastic = true
    else
        maxits = 30
        Δλ     = 0.0
        # Δλ     = abs(τtr)/ks
        τy     = 0.0
        tol    = 1e-4


        for i in 1:maxits
            # τ  = τtr - Δλ*ks*sign(τtr)
            τy = Tau(mat, ipd.sy+Δλ)
            # τ  = τy*sign(τtr)
            dτydsy = deriv(mat, ipd, ipd.sy+Δλ)

            # @s τ
            # @s τy
            # @s dτydsy

            # f     = abs(τ) - τy
            f = abs(τtr) - Δλ*ks - τy
            # dfdΔλ = -sign(abs(τtr) - Δλ*ks)*ks - dτydsy
            dfdΔλ = -ks - dτydsy
            Δλ = Δλ - f/dfdΔλ
            # Δλ = (abs(τtr)-τy)/mat.ks
            # abs(f) < tol && @s i
            abs(f) < tol && break

            # @s f
            # @s dfdΔλ

            # @s Δλ
            # if Δλ<0
                # Δλ = (abs(τtr)-ipd.τy)/mat.ks
                # break
            # end

            if i==maxits || isnan(Δλ)  
                if Δλ<0.0
                    @s Δλ
                    Δλ = (abs(τtr)-ipd.τy)/mat.ks
                    break
                end

                               
                @s "Errorrrrr"
                @s dτydsy

                # error()
                # @s f
                # @s i
                return ipd.σ, failure("CEBJoint1D: Could not find Δλ")
                # break
            end

        end

        if Δλ<0.0
            Δλ = (abs(τtr)-ipd.τy)/mat.ks
            @s Δλ
        end

        ipd.sy += Δλ
        ipd.τy  = τy
        τ  = ipd.τy*sign(τtr)
        Δτ = τ - τini
        ipd.elastic = false
    end

    # calculate Δσ
    Δτ = τ - τini
    Δσ = kn*Δu
    Δσ[1] = Δτ

    # update u and σ
    ipd.u .+= Δu
    ipd.σ .+= Δσ

    return Δσ, success()
end


function stress_update(mat::CebRSJoint, ipd::CebRSJointState, Δu::Vect)
    ks = mat.ks
    kn = mat.kn
    Δs = Δu[1]      # relative displacement
    τini = ipd.σ[1] # initial shear stress
    τtr  = τini + ks*Δs # elastic trial

    ftr  = yield_func(mat, ipd, τtr)

    if ftr<0.0
        τ = τtr
        ipd.elastic = true
    else
        dτydsy = deriv(mat, ipd, ipd.sy)
        # @s ks
        # @s dτydsy
        # Δsy     = (abs(τtr)-ipd.τy)/(ks+dτydsy)
        # Δsy     = (abs(τtr)-ipd.τy)/abs(ks+dτydsy)
        Δsy     = (abs(τtr)-ipd.τy)/ks
        ipd.sy += Δsy
        ipd.τy  = Tau(mat, ipd.sy)
        τ  = ipd.τy*sign(τtr)
        Δτ = τ - τini
        ipd.elastic = false
    end

    # calculate Δσ
    Δτ = τ - τini
    Δσ = kn*Δu
    Δσ[1] = Δτ

    # update u and σ
    ipd.u .+= Δu
    ipd.σ .+= Δσ

    return Δσ, success()
end

function ip_state_vals(mat::CebRSJoint, ipd::CebRSJointState)
    return OrderedDict(
      :ur   => ipd.u[1] ,
      :tau  => ipd.σ[1] ,
      )
end

