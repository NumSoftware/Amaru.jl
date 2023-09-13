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
    ν::Float64

    function ElasticBeam(; params...)
        names = (E="Young modulus", nu="Poisson ratio")
        required = (:E, :nu)
        @checkmissing params required names

        params  = values(params)
        E       = params.E
        nu      = params.nu

        @check E>=0.0
        @check 0<=nu<0.5
        return new(E, nu)
    end
end


# Type of corresponding state structure
compat_state_type(::Type{ElasticBeam}) = ElasticBeamState

# Element types that work with this material
compat_elem_types(::Type{ElasticBeam}) = (MechBeam,)


function calcD(mat::ElasticBeam, state::ElasticBeamState)
    E = mat.E
    ν = mat.ν
    c = E/(1-ν^2)
    g = E/(1+ν)

    if state.env.ndim==2
        return [ c      0.0
                 0.0  5/6*g ]
    else
        return [ c   0.0    0.0  
                0.0  5/6*g  0.0  
                0.0  0.0    5/6*g ]
    end
end


function update_state!(mat::ElasticBeam, state::ElasticBeamState, dε::Array{Float64,1})
    D = calcD(mat, state)
    dσ = D*dε
    state.ε += dε
    state.σ += dσ
    return dσ, success()
end


function ip_state_vals(mat::ElasticBeam, state::ElasticBeamState)
    vals = OrderedDict{Symbol,Float64}(
      :sX  => state.σ[1],
      :eX  => state.ε[1],
      :sXY => state.σ[2]/SR2
    )
    if state.env.ndim==3
        vals[:sXZ] = state.σ[2]/SR2
        vals[:sXY] = state.σ[3]/SR2
    end
    return vals
end
