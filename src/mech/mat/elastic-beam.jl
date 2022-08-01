# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticBeam

mutable struct ElasticBeamState<:IpState
    "Environment information"
    env::ModelEnv
    "Stress tensor"
    σ::Array{Float64,1} # in the local system
    "Strain tensor"
    ε::Array{Float64,1} # in the local system

    function ElasticBeamState(env::ModelEnv=ModelEnv())
        this = new(env)
        nstr = env.ndim==2 ? 2 : 3
        this.σ = zeros(nstr)
        this.ε = zeros(nstr)
        return this
    end
end

mutable struct ElasticBeam<:Material
    E::Float64
    nu::Float64
    thz::Float64
    thy::Float64
    A::Float64
    ρ::Float64

    function ElasticBeam(prms::Dict{Symbol,Float64})
        return  ElasticBeam(;prms...)
    end

    function ElasticBeam(;E=NaN, nu=0.0, thz=NaN, thy=NaN, rho=0.0)
        @check E>0.0
        @check 0.0<=nu<0.5
        @check thz>0.0
        @check thy>0.0
        @check rho>=0.0
        this = new(E, nu, thy, thz, thy*thz, rho)
        return this
    end
end

matching_elem_type(::ElasticBeam) = MechBeam

# Type of corresponding state structure
ip_state_type(mat::ElasticBeam) = ElasticBeamState


function calcD(mat::ElasticBeam, state::ElasticBeamState)
    E = mat.E
    ν = mat.nu
    c = E/(1.0-ν^2)
    g = E/(1+ν)

    if state.env.ndim==2
        return [ c      0.0  
                 0.0  5/6*g ]
    else
        return [ c    0.0    0.0  
                0.0  5/6*g  0.0  
                0.0  0.0    5/6*g ]
    end
end


function stress_update(mat::ElasticBeam, state::ElasticBeamState, dε::Array{Float64,1})
    D = calcD(mat, state)
    dσ = D*dε
    state.ε += dε
    state.σ += dσ
    return dσ, success()
end


function ip_state_vals(mat::ElasticBeam, state::ElasticBeamState)
    vals =  OrderedDict(
      "sx'"   => state.σ[1],
      "sx'y'" => state.σ[2],
      "ex'"   => state.ε[1],
      "A"     => mat.A )
    if state.env.ndim==3
        vals["sx'z'"] = state.σ[2]
        vals["sx'y'"] = state.σ[3]
    end
    return vals
end
