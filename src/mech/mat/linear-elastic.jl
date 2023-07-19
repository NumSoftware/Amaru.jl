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
    ν::Float64

    @doc """
        $(SIGNATURES)

    Creates an `LinearElastic` material type

    # Arguments
    - `E`: Young modulus
    - `nu`: Poisson ratio
    """
    function LinearElastic(; params...)
        names = (E="Young modulus", nu="Poisson ratio")
        required = (:E, :nu)
        @checkmissing params required names

        params  = values(params)
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


# Type of corresponding state structure
ip_state_type(::Type{LinearElastic}) = ElasticSolidState

# Element types that work with this material
matching_elem_types(::Type{LinearElastic}) = (MechSolid, MechShell)


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
    elseif stressmodel=="shell"
        c = E/(1.0-ν^2)
        g = E/(1+ν)
        return @Mat6x6 [
            c    c*ν   0.0  0.0    0.0    0.0
            c*ν  c     0.0  0.0    0.0    0.0
            0.0  0.0   0.0  0.0    0.0    0.0
            0.0  0.0   0.0  5/6*g  0.0    0.0
            0.0  0.0   0.0  0.0    5/6*g  0.0
            0.0  0.0   0.0  0.0    0.0    g ]
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


function calcD(mat::LinearElastic, state::ElasticSolidState, stressmodel::String=state.env.ana.stressmodel)
    return calcDe(mat.E, mat.ν, stressmodel)
end


function update_state!(mat::LinearElastic, state::ElasticSolidState, dε::AbstractArray, stressmodel::String=state.env.ana.stressmodel)
    De = calcDe(mat.E, mat.ν, stressmodel)
    dσ = De*dε
    state.ε += dε
    state.σ += dσ
    return dσ, success()
end


function ip_state_vals(mat::LinearElastic, state::ElasticSolidState, stressmodel::String=state.env.ana.stressmodel)
    return stress_strain_dict(state.σ, state.ε, stressmodel)
end