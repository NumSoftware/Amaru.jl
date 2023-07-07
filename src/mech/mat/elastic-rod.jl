# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticRod

"""
    ElasticRod

A type for linear elastic materials in rods.

# Fields

$(TYPEDFIELDS)
"""
mutable struct ElasticRod<:MatParams
    "Young Modulus"
    E::Float64
    # "Section area"
    # A::Float64
    # "Density"
    # ρ::Float64

    function ElasticRod(prms::Dict{Symbol,Float64})
        return  ElasticRod(;prms...)
    end

    @doc """
        $(SIGNATURES)

    Creates an `ElasticRod` material type

    # Arguments
    - `E`: Young modulus
    - `A`: Section area
    - `dm`: Diameter (only if `A` is not provided)
    - `rho`: Density
    """
    function ElasticRod(;E::Number=NaN)
        @check E>0.0
        return new(E)
        # @check A>0.0 || dm>0.0
        # @check rho>=0.0
        # dm>0 && (A=π*dm^2/4)
        # return new(E, A, rho)
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

matching_elem_type(::ElasticRod) = MechRodElem
matching_elem_type_if_embedded(::ElasticRod) = MechEmbRodElem

# Type of corresponding state structure
ip_state_type(matparams::ElasticRod) = ElasticRodState


function update_state(matparams::ElasticRod, state::ElasticRodState, Δε::Float64)
    Δσ = matparams.E*Δε
    state.ε += Δε
    state.σ += Δσ
    return Δσ, success()
end


function ip_state_vals(matparams::ElasticRod, state::ElasticRodState)
    return OrderedDict(
      :sa => state.σ,
      :ea => state.ε,
    #   :fa => state.σ*matparams.A,
    #   :A  => matparams.A 
      )
end


function calcD(matparams::ElasticRod, ips::ElasticRodState)
    return matparams.E
end

