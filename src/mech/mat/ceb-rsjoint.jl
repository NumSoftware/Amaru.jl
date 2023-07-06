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


mutable struct CebRSJoint<:MatParams
    τmax:: Float64
    τres:: Float64
    s1  :: Float64
    s2  :: Float64
    s3  :: Float64
    α   :: Float64
    β   :: Float64
    ks  :: Float64
    kn  :: Float64

    function CebRSJoint(prms::Dict{Symbol,Float64})
        return  CebRSJoint(;prms...)
    end

    function CebRSJoint(;taumax=NaN, taures=NaN, TauM=NaN, TauR=NaN, s1=NaN, s2=NaN, s3=NaN, alpha=NaN, beta=NaN, kn=NaN, ks=NaN)
        
        @check s1>0
        @check s2>s1
        @check s3>s2

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

        this = new(taumax, taures, s1, s2, s3, alpha, beta, ks, kn)
        return this
    end
end

# Returns the element type that works with this material
matching_elem_type(::CebRSJoint) = MechRSJointElem

# Creates a new instance of Ip data
ip_state_type(::MechRSJointElem, ::CebRSJoint) = CebRSJointState

CEBJoint1D = CebRSJoint #! deprecated
export CEBJoint1D


function Tau(matparams::CebRSJoint, sy::Float64)
    if sy<matparams.s1
        return matparams.τmax*(sy/matparams.s1)^matparams.α
    elseif sy<matparams.s2
        return matparams.τmax
    elseif sy<matparams.s3
        return matparams.τmax - (matparams.τmax-matparams.τres)*((sy-matparams.s2)/(matparams.s3-matparams.s2))^matparams.β
    else
        return matparams.τres  + 1*(sy-matparams.s3)
    end
end


function deriv(matparams::CebRSJoint, state::CebRSJointState, sy::Float64)
    if sy==0.0
        s1_factor = 0.01
        sy = s1_factor*matparams.s1   # to avoid undefined derivative
    end

    if sy<=matparams.s1
        return matparams.τmax/matparams.s1*(sy/matparams.s1)^(matparams.α-1)
    elseif sy<matparams.s2
        #return matparams.ks/10000
        return 1.0
    elseif sy<matparams.s3
        return -(matparams.τmax-matparams.τres)/(matparams.s3-matparams.s2)*((sy-matparams.s2)/(matparams.s3-matparams.s2))^(matparams.β-1)
    else
        return 1.0
        # return matparams.ks/10000
    end
end


function calcD(matparams::CebRSJoint, state::CebRSJointState)
    ndim = state.env.ndim
    ks = matparams.ks

    if !state.elastic
        dτydsy = deriv(matparams, state, state.sy)
        # ks = ks*dτydsy/(ks+dτydsy)
        ks = dτydsy
        # @s ks
    end

    kn = matparams.kn
    if state.env.ndim==2
        return [  ks  0.0
                 0.0   kn ]
    else
        return  [  ks  0.0  0.0
                  0.0   kn  0.0
                  0.0  0.0   kn ]
    end
end


function yield_func(matparams::CebRSJoint, state::CebRSJointState, τ::Float64)
    return abs(τ) - state.τy
end


function stress_update_n(matparams::CebRSJoint, state::CebRSJointState, Δu::Vect)
    ks = matparams.ks
    kn = matparams.kn
    Δs = Δu[1]      # relative displacement
    τini = state.σ[1] # initial shear stress
    τtr  = τini + ks*Δs # elastic trial

    ftr  = yield_func(matparams, state, τtr)

    if ftr<0.0
        τ = τtr
        state.elastic = true
    else
        maxits = 30
        Δλ     = 0.0
        # Δλ     = abs(τtr)/ks
        τy     = 0.0
        tol    = 1e-4


        for i in 1:maxits
            # τ  = τtr - Δλ*ks*sign(τtr)
            τy = Tau(matparams, state.sy+Δλ)
            # τ  = τy*sign(τtr)
            dτydsy = deriv(matparams, state, state.sy+Δλ)

            # @s τ
            # @s τy
            # @s dτydsy

            # f     = abs(τ) - τy
            f = abs(τtr) - Δλ*ks - τy
            # dfdΔλ = -sign(abs(τtr) - Δλ*ks)*ks - dτydsy
            dfdΔλ = -ks - dτydsy
            Δλ = Δλ - f/dfdΔλ
            # Δλ = (abs(τtr)-τy)/matparams.ks
            # abs(f) < tol && @s i
            abs(f) < tol && break

            # @s f
            # @s dfdΔλ

            # @s Δλ
            # if Δλ<0
                # Δλ = (abs(τtr)-state.τy)/matparams.ks
                # break
            # end

            if i==maxits || isnan(Δλ)  
                if Δλ<0.0
                    @s Δλ
                    Δλ = (abs(τtr)-state.τy)/matparams.ks
                    break
                end

                               
                @s "Errorrrrr"
                @s dτydsy

                # error()
                # @s f
                # @s i
                return state.σ, failure("CEBJoint1D: Could not find Δλ")
                # break
            end

        end

        if Δλ<0.0
            Δλ = (abs(τtr)-state.τy)/matparams.ks
            @s Δλ
        end

        state.sy += Δλ
        state.τy  = τy
        τ  = state.τy*sign(τtr)
        Δτ = τ - τini
        state.elastic = false
    end

    # calculate Δσ
    Δτ = τ - τini
    Δσ = kn*Δu
    Δσ[1] = Δτ

    # update u and σ
    state.u .+= Δu
    state.σ .+= Δσ

    return Δσ, success()
end


function update_state(matparams::CebRSJoint, state::CebRSJointState, Δu::Vect)
    ks = matparams.ks
    kn = matparams.kn
    Δs = Δu[1]      # relative displacement
    τini = state.σ[1] # initial shear stress
    τtr  = τini + ks*Δs # elastic trial

    ftr  = yield_func(matparams, state, τtr)

    if ftr<0.0
        τ = τtr
        state.elastic = true
    else
        dτydsy = deriv(matparams, state, state.sy)
        # @s ks
        # @s dτydsy
        # Δsy     = (abs(τtr)-state.τy)/(ks+dτydsy)
        # Δsy     = (abs(τtr)-state.τy)/abs(ks+dτydsy)
        Δsy     = (abs(τtr)-state.τy)/ks
        state.sy += Δsy
        state.τy  = Tau(matparams, state.sy)
        τ  = state.τy*sign(τtr)
        Δτ = τ - τini
        state.elastic = false
    end

    # calculate Δσ
    Δτ = τ - τini
    Δσ = kn*Δu
    Δσ[1] = Δτ

    # update u and σ
    state.u .+= Δu
    state.σ .+= Δσ

    return Δσ, success()
end

function ip_state_vals(matparams::CebRSJoint, state::CebRSJointState)
    return OrderedDict(
      :ur   => state.u[1] ,
      :tau  => state.σ[1] ,
      )
end

