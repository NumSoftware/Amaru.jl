# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LinearElastic

"""
    LinearElastic

A type for linear elastic materials.

# Fields

$(TYPEDFIELDS)
"""
mutable struct LinearElastic<:MatParams
    "Young Modulus"
    E ::Float64
    "Poisson ratio"
    nu::Float64
    "Density"
    ρ::Float64

    function LinearElastic(prms::Dict{Symbol,Float64})
        return  LinearElastic(;prms...)
    end

    @doc """
        $(SIGNATURES)

    Creates an `LinearElastic` material type

    # Arguments
    - `E`: Young modulus
    - `nu`: Poisson ratio
    - `rho`: Density
    """
    function LinearElastic(;E=1.0, nu=0.0, rho=0.0)
        E<=0.0       && error("Invalid value for E: $E")
        !(0<=nu<0.5) && error("Invalid value for nu: $nu")
        rho<0.0      && error("Invalid value for rho: $rho")
        this = new(E, nu, rho)
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
    σ::Array{Float64,1}
    "Strain tensor"
    ε::Array{Float64,1}

    function ElasticSolidState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ = zeros(6)
        this.ε = zeros(6)
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


mutable struct ElasticShellState<:IpState
    "Environment information"
    env::ModelEnv
    "Stress tensor"
    σ::Array{Float64,1} # in the local system
    "Strain tensor"
    ε::Array{Float64,1} # in the local system

    function ElasticShellState(env::ModelEnv=ModelEnv())
        this = new(env)
        this.σ = zeros(6)
        this.ε = zeros(6)
        return this
    end
end


# Returns the default element type that works with this material model



# Type of corresponding state structure
ip_state_type(::MechSolidElem, ::LinearElastic) = ElasticSolidState
ip_state_type(::MechRodElem, ::LinearElastic)  = ElasticRodState
ip_state_type(::MechBeamElem, ::LinearElastic) = ElasticBeamState
ip_state_type(::MechShellElem, ::LinearElastic) = ElasticShellState


function calcDe(E::Number, ν::Number, stressmodel::String)
    if stressmodel=="plane-stress"
        c = E/(1.0-ν^2)
        return [
            c    c*ν   0.0  0.0  0.0  0.0
            c*ν  c     0.0  0.0  0.0  0.0
            0.0  0.0   0.0  0.0  0.0  0.0
            0.0  0.0   0.0  0.0  0.0  0.0
            0.0  0.0   0.0  0.0  0.0  0.0
            0.0  0.0   0.0  0.0  0.0  c*(1.0-ν) ]
        ezz = -ν/E*(sxx+syy)
    else
        c = E/((1.0+ν)*(1.0-2.0*ν))
        return [
            c*(1-ν) c*ν     c*ν     0.0         0.0         0.0
            c*ν     c*(1-ν) c*ν     0.0         0.0         0.0
            c*ν     c*ν     c*(1-ν) 0.0         0.0         0.0
            0.0     0.0     0.0     c*(1-2*ν)   0.0         0.0
            0.0     0.0     0.0     0.0         c*(1-2*ν)   0.0
            0.0     0.0     0.0     0.0         0.0         c*(1-2*ν) ]
    end
end


function calcD(matparams::LinearElastic, state::ElasticSolidState)
    return calcDe(matparams.E, matparams.nu, state.env.anaprops.stressmodel)
end


function update_state(matparams::LinearElastic, state::ElasticSolidState, dε::Array{Float64,1})
    De = calcDe(matparams.E, matparams.nu, state.env.anaprops.stressmodel)
    dσ = De*dε
    state.ε += dε
    state.σ += dσ
    return dσ, success()
end


function ip_state_vals(matparams::LinearElastic, state::ElasticSolidState)
    return stress_strain_dict(state.σ, state.ε, state.env.anaprops.stressmodel)
end

# Rod


function update_state(matparams::LinearElastic, state::ElasticRodState, Δε::Float64)
    Δσ = matparams.E*Δε
    state.ε += Δε
    state.σ += Δσ
    return Δσ, success()
end


function ip_state_vals(matparams::LinearElastic, state::ElasticRodState)
    return OrderedDict(
      :sa => state.σ,
      :ea => state.ε,
      )
end


function calcD(matparams::LinearElastic, ips::ElasticRodState)
    return matparams.E
end

# Beam


function calcD(matparams::LinearElastic, state::ElasticBeamState)
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


function update_state(matparams::LinearElastic, state::ElasticBeamState, dε::Array{Float64,1})
    D = calcD(matparams, state)
    dσ = D*dε
    state.ε += dε
    state.σ += dσ
    return dσ, success()
end


function ip_state_vals(matparams::LinearElastic, state::ElasticBeamState)
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

# Shell


function calcD(matparams::LinearElastic, state::ElasticShellState)
    E = matparams.E
    ν = matparams.nu
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


function update_state(matparams::LinearElastic, state::ElasticShellState, dε::Array{Float64,1})
    D = calcD(matparams, state)
    dσ = D*dε
    state.ε += dε
    state.σ += dσ
    return dσ, success()
end


function ip_state_vals(matparams::LinearElastic, state::ElasticShellState)
    return stress_strain_dict(state.σ, state.ε, state.env.anaprops.stressmodel)
end