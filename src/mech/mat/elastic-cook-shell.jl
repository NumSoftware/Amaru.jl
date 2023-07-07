# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticCookShell

mutable struct ElasticCookShellState<:IpState
    "Environment information"
    env::ModelEnv
    "Stress tensor"
    σ::Array{Float64,1}
    "Strain tensor"
    ε::Array{Float64,1}

    function ElasticCookShellState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ = zeros(6)
        this.ε = zeros(6)
        return this
    end
end

mutable struct ElasticCookShell<:MatParams
    E::Float64
    nu::Float64
    th::Float64

    function ElasticCookShell(prms::Dict{Symbol,Float64})
        return  ElasticCookShell(;prms...)
    end

    function ElasticCookShell(;E=NaN, nu=NaN, thickness=0.1)
        E>0.0 || error("Invalid value for E: $E")
        (0<=nu<0.5) || error("Invalid value for nu: $nu")
        thickness >0.0 || error("Invalid value for thickness: $thickness")
        this = new(E, nu, thickness)
        return this
    end
end

matching_elem_type(::ElasticCookShell) = CookShellElem

# Type of corresponding state structure
ip_state_type(matparams::ElasticCookShell) = ElasticCookShellState


function calcD(matparams::ElasticCookShell, state::ElasticCookShellState)
    E = matparams.E
    ν = matparams.nu
    c = E/(1.0-ν^2)
    # g = c*(1.0-ν)
    g = E/(1+ν)
    return [
        c    c*ν   0.0  0.0    0.0    0.0
        c*ν  c     0.0  0.0    0.0    0.0
        0.0  0.0   0.0  0.0    0.0    0.0
        0.0  0.0   0.0  5/6*g  0.0    0.0
        0.0  0.0   0.0  0.0    5/6*g  0.0
        0.0  0.0   0.0  0.0    0.0    g ]
    # ezz = -ν/E*(sxx+syy)
end


function update_state(matparams::ElasticCookShell, state::ElasticCookShellState, dε::Array{Float64,1})
    D = calcD(matparams, state)
    dσ = D*dε
    state.ε += dε
    state.σ += dσ
    return dσ, success()
end


function ip_state_vals(matparams::ElasticCookShell, state::ElasticCookShellState)
    return stress_strain_dict(state.σ, state.ε, state.env.anaprops.stressmodel)
end