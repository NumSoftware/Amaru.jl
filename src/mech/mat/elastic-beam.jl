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

mutable struct ElasticBeam<:MatParams
    E::Float64
    nu::Float64

    function ElasticBeam(prms::Dict{Symbol,Float64})
        return  ElasticBeam(;prms...)
    end

    function ElasticBeam(;E=NaN, nu=0.0)
        @check E>0.0
        @check 0.0<=nu<0.5
        this = new(E, nu)
        return this
    end
end

matching_elem_type(::ElasticBeam) = MechBeamElem

# Type of corresponding state structure
ip_state_type(matparams::ElasticBeam) = ElasticBeamState


function calcD(matparams::ElasticBeam, state::ElasticBeamState)
    E = matparams.E
    ν = matparams.nu
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


function update_state(matparams::ElasticBeam, state::ElasticBeamState, dε::Array{Float64,1})
    D = calcD(matparams, state)
    dσ = D*dε
    state.ε += dε
    state.σ += dσ
    return dσ, success()
end


function ip_state_vals(matparams::ElasticBeam, state::ElasticBeamState)
    vals =  OrderedDict(
      "sx'"   => state.σ[1],
      "ex'"   => state.ε[1],
      "sx'y'" => state.σ[2]/SR2)
    if state.env.ndim==3
        vals["sx'z'"] = state.σ[2]/SR2
        vals["sx'y'"] = state.σ[3]/SR2
    end
    return vals
end
