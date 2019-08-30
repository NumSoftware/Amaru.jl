# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export CEBJoint1D

mutable struct CEBJoint1DIpState<:IpState
    env::ModelEnv
    σ  ::Array{Float64,1}
    u  ::Array{Float64,1}
    τy ::Float64      # max stress
    sy ::Float64      # accumulated relative displacement
    elastic::Bool
    function CEBJoint1DIpState(env::ModelEnv=ModelEnv())
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

mutable struct CEBJoint1D<:Material
    τmax:: Float64
    τres:: Float64
    s1  :: Float64
    s2  :: Float64
    s3  :: Float64
    α   :: Float64
    β   :: Float64
    ks  :: Float64
    kn  :: Float64
    h   :: Float64

    function CEBJoint1D(prms::Dict{Symbol,Float64})
        return  CEBJoint1D(;prms...)
    end

    function CEBJoint1D(;TauM=NaN, TauR=NaN, s1=NaN, s2=NaN, s3=NaN, alpha=NaN, beta=NaN, kn=NaN, ks=NaN, h=NaN, A=NaN, dm=NaN)
        @assert s1>0
        @assert s2>s1
        @assert s3>s2
        @assert ks>0

        # Estimate the perimeter h
        if isnan(h)
            if A>0
                h = 2.0*√(A*pi)
            else
                @assert dm>0
                h = pi*dm
            end
        end
        @assert h>0

        # Estimate TauM if not provided
        if isnan(TauM)
            TauM = ks*s1
        end
        @assert TauM>TauR

        # Define alpha if not provided
        if isnan(alpha); alpha = 1.0 end
        @assert 0.0<=alpha<=1.0

        # Define beta if not provided
        if isnan(beta); beta= 1.0 end
        @assert beta>=0.0

        # Estimate kn if not provided
        if isnan(kn)
            kn = ks
        end
        @assert kn>0

        this = new(TauM, TauR, s1, s2, s3, alpha, beta, ks, kn, h)
        return this
    end
end

# Returns the element type that works with this material
matching_elem_type(::CEBJoint1D) = MechJoint1D

# Creates a new instance of Ip data
ip_state_type(mat::CEBJoint1D) = CEBJoint1DIpState


function Tau(mat::CEBJoint1D, sy::Float64)
    if sy<mat.s1
        return mat.τmax*(sy/mat.s1)^mat.α
    elseif sy<mat.s2
        return mat.τmax
    elseif sy<mat.s3
        return mat.τmax - (mat.τmax-mat.τres)*((sy-mat.s2)/(mat.s3-mat.s2))^mat.β
    else
        return mat.τres
    end
end

function deriv(mat::CEBJoint1D, ipd::CEBJoint1DIpState, sy::Float64)
    if sy==0.0
        S1_FACTOR = 0.01
        sy = S1_FACTOR*mat.s1   # to avoid undefined derivative
    end

    if sy<=mat.s1
        return mat.τmax/mat.s1*(sy/mat.s1)^(mat.α-1)
    elseif sy<mat.s2
        return mat.ks/1000
        #return 0.0
    elseif sy<mat.s3
        return -(mat.τmax-mat.τres)/(mat.s3-mat.s2)*((sy-mat.s2)/(mat.s3-mat.s2))^(mat.β-1)
    else
        #return 0.0
        return mat.ks/1000
    end
end

function calcD(mat::CEBJoint1D, ipd::CEBJoint1DIpState)
    ndim = ipd.env.ndim
    if ipd.elastic
        ks = mat.ks
    else
        ks = deriv(mat, ipd, ipd.sy)
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

function yield_func(mat::CEBJoint1D, ipd::CEBJoint1DIpState, τ::Float64)
    return abs(τ) - ipd.τy
end

function stress_update(mat::CEBJoint1D, ipd::CEBJoint1DIpState, Δu::Vect)
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
        Δsy     = (abs(τtr)-ipd.τy)/mat.ks
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

    return Δσ
end

function ip_state_vals(mat::CEBJoint1D, ipd::CEBJoint1DIpState)
    return OrderedDict(
      :ur   => ipd.u[1] ,
      :tau  => ipd.σ[1] ,
      )
end

