# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LinearElastic

"""
    LinearElastic

A type for linear elastic materials.

# Fields

$(TYPEDFIELDS)
"""
mutable struct LinearElastic<:Material
    "Young Modulus"
    E ::Float64
    "Poisson ratio"
    nu::Float64

    @doc """
        $(SIGNATURES)

    Creates an `LinearElastic` material type

    # Arguments
    - `E`: Young modulus
    - `nu`: Poisson ratio
    """
    function LinearElastic(;params...)
        names = (E="Young modulus", nu="Poisson ratio")
        required = (:E,)
        # @show params
        @checkmissing params required names

        default = (nu=0.0,)
        params  = merge(default, params)
        E       = params.E
        nu      = params.nu

        @check E>=0.0
        @check 0<=nu<0.5
        this = new(E, nu)
        return this
    end
end


"""
    ElasticSolidState

A type for the state data of a `LinearElastic` type.

# Fields

$(TYPEDFIELDS)
"""
mutable struct ElasticSolidState<:IpState
    "Environment information"
    env::ModelEnv
    "Stress tensor"
    σ::Vec6
    "Strain tensor"
    ε::Vec6

    function ElasticSolidState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ = zeros(Vec6)
        this.ε = zeros(Vec6)
        return this
    end
end


"""
    ElasticRodState

A type for the state data of a `ElasticRod` type.

# Fields

$(TYPEDFIELDS)
"""
mutable struct ElasticRodState<:IpState
    "environment information"
    env::ModelEnv
    "Axial stress"
    σ::Float64
    "Axial strain"
    ε::Float64
    function ElasticRodState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ = 0.0
        this.ε = 0.0
        return this
    end
end


mutable struct ElasticBeamState<:IpState
    "Environment information"
    env::ModelEnv
    "Stress tensor"
    σ::Vect # in the local system
    "Strain tensor"
    ε::Vect # in the local system

    function ElasticBeamState(env::ModelEnv=ModelEnv())
        this = new(env)
        nstr = env.ndim==2 ? 2 : 3
        this.σ = zeros(nstr)
        this.ε = zeros(nstr)
        return this
    end
end


mutable struct ElasticShellState<:IpState
    "Environment information"
    env::ModelEnv
    "Stress tensor"
    σ::Vec6 # in the local system
    "Strain tensor"
    ε::Vec6 # in the local system

    function ElasticShellState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ = zeros(6)
        this.ε = zeros(6)
        return this
    end
end


# Type of corresponding state structure
ip_state_type(::MechSolid, ::LinearElastic) = ElasticSolidState
ip_state_type(::MechRod, ::LinearElastic)  = ElasticRodState
ip_state_type(::MechBeam, ::LinearElastic) = ElasticBeamState
ip_state_type(::MechShell, ::LinearElastic) = ElasticShellState


function calcDe(E::Number, ν::Number, stressmodel::String)
    if stressmodel=="plane-stress"
        c = E/(1.0-ν^2)
        return @Mat6x6 [
            c    c*ν   0.0  0.0  0.0  0.0
            c*ν  c     0.0  0.0  0.0  0.0
            0.0  0.0   0.0  0.0  0.0  0.0
            0.0  0.0   0.0  0.0  0.0  0.0
            0.0  0.0   0.0  0.0  0.0  0.0
            0.0  0.0   0.0  0.0  0.0  c*(1.0-ν) ]
        ezz = -ν/E*(sxx+syy)
    else
        c = E/((1.0+ν)*(1.0-2.0*ν))
        return @Mat6x6 [
            c*(1-ν) c*ν     c*ν     0.0         0.0         0.0
            c*ν     c*(1-ν) c*ν     0.0         0.0         0.0
            c*ν     c*ν     c*(1-ν) 0.0         0.0         0.0
            0.0     0.0     0.0     c*(1-2*ν)   0.0         0.0
            0.0     0.0     0.0     0.0         c*(1-2*ν)   0.0
            0.0     0.0     0.0     0.0         0.0         c*(1-2*ν) ]
    end
end


function calcD(mat::LinearElastic, state::ElasticSolidState)
    return calcDe(mat.E, mat.nu, state.env.ana.stressmodel)
end


function update_state(mat::LinearElastic, state::ElasticSolidState, dε::AbstractArray)
    De = calcDe(mat.E, mat.nu, state.env.ana.stressmodel)
    dσ = De*dε
    state.ε += dε
    state.σ += dσ
    return dσ, success()
end


function ip_state_vals(mat::LinearElastic, state::ElasticSolidState)
    return stress_strain_dict(state.σ, state.ε, state.env.ana.stressmodel)
end


## Rod


function update_state(mat::LinearElastic, state::ElasticRodState, Δε::Float64)
    Δσ = mat.E*Δε
    state.ε += Δε
    state.σ += Δσ
    return Δσ, success()
end


function ip_state_vals(mat::LinearElastic, state::ElasticRodState)
    return OrderedDict(
      :sa => state.σ,
      :ea => state.ε,
      )
end


function calcD(mat::LinearElastic, ips::ElasticRodState)
    return mat.E
end


## Beam


function calcD(mat::LinearElastic, state::ElasticBeamState)
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


function update_state(mat::LinearElastic, state::ElasticBeamState, dε::AbstractArray)
    D = calcD(mat, state)
    dσ = D*dε
    state.ε += dε
    state.σ += dσ
    return dσ, success()
end


function ip_state_vals(mat::LinearElastic, state::ElasticBeamState)
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


## Shell


function calcD(mat::LinearElastic, state::ElasticShellState)
    E = mat.E
    ν = mat.nu
    c = E/(1.0-ν^2)
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


function update_state(mat::LinearElastic, state::ElasticShellState, dε::AbstractArray)
    D = calcD(mat, state)
    dσ = D*dε
    state.ε += dε
    state.σ += dσ
    return dσ, success()
end


function ip_state_vals(mat::LinearElastic, state::ElasticShellState)
    return stress_strain_dict(state.σ, state.ε, state.env.ana.stressmodel)
end