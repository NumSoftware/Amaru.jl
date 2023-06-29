# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticSolid

"""
    ElasticSolid

A type for linear elastic materials.

# Fields

$(TYPEDFIELDS)
"""
mutable struct ElasticSolid<:Material
    "Young Modulus"
    E ::Float64
    "Poisson ratio"
    nu::Float64
    "Density"
    ρ::Float64

    function ElasticSolid(prms::Dict{Symbol,Float64})
        return  ElasticSolid(;prms...)
    end

    @doc """
        $(SIGNATURES)

    Creates an `ElasticSolid` material type

    # Arguments
    - `E`: Young modulus
    - `nu`: Poisson ratio
    - `rho`: Density
    """
    function ElasticSolid(;E=1.0, nu=0.0, rho=0.0)
        E<=0.0       && error("Invalid value for E: $E")
        !(0<=nu<0.5) && error("Invalid value for nu: $nu")
        rho<0.0      && error("Invalid value for rho: $rho")
        this = new(E, nu, rho)
        return this
    end
end

"""
    ElasticSolidState

A type for the state data of a `ElasticSolid` type.

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


# Returns the element type that works with this material model
function matching_elem_type(::ElasticSolid, shape::CellShape, ndim::Int) 
    if shape.ndim==2 && ndim==3
        return MechShell
    else
        return MechSolid
    end
end

# Type of corresponding state structure
ip_state_type(::ElasticSolid) = ElasticSolidState


function calcDe(E::Number, ν::Number, modeltype::String)
    if modeltype=="plane-stress"
        c = E/(1.0-ν^2)
        return [
            c    c*ν   0.0  0.0  0.0  0.0
            c*ν  c     0.0  0.0  0.0  0.0
            0.0  0.0   0.0  0.0  0.0  0.0
            0.0  0.0   0.0  0.0  0.0  0.0
            0.0  0.0   0.0  0.0  0.0  0.0
            0.0  0.0   0.0  0.0  0.0  c*(1.0-ν) ]
        # ezz = -ν/E*(sxx+syy)
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

function calcD(mat::ElasticSolid, state::ElasticSolidState)
    return calcDe(mat.E, mat.nu, state.env.modeltype)
end

function stress_update(mat::ElasticSolid, state::ElasticSolidState, dε::Array{Float64,1})
    De = calcDe(mat.E, mat.nu, state.env.modeltype)
    dσ = De*dε
    state.ε += dε
    state.σ += dσ
    return dσ, success()
end

function ip_state_vals(mat::ElasticSolid, state::ElasticSolidState)
    return stress_strain_dict(state.σ, state.ε, state.env.modeltype)
end
