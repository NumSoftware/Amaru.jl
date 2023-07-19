# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export ElasticRod

"""
    ElasticRod

A type for linear elastic materials in rods.

# Fields

$(TYPEDFIELDS)
"""
mutable struct ElasticRod<:Material
    "Young Modulus"
    E::Float64

    @doc """
        $(SIGNATURES)

    Creates an `ElasticRod` material type

    # Arguments
    - `E`: Young modulus
    - `rho`: Density
    """
    function ElasticRod(; params...)
        names = (E="Young modulus")
        required = (:E,)
        @checkmissing params required names

        params  = values(params)
        E       = params.E

        @check E>=0.0
        return new(E)
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



# Type of corresponding state structure
ip_state_type(::Type{ElasticRod}) = ElasticRodState

# Element types that work with this material
matching_elem_types(::Type{ElasticRod}) = (MechRod,)
matching_elem_type_if_embedded(::Type{ElasticRod}) = MechEmbRod


function update_state!(mat::ElasticRod, state::ElasticRodState, Δε::Float64)
    Δσ = mat.E*Δε
    state.ε += Δε
    state.σ += Δσ
    return Δσ, success()
end


function ip_state_vals(mat::ElasticRod, state::ElasticRodState)
    return OrderedDict(
      :sa => state.σ,
      :ea => state.ε,
    #   :fa => state.σ*mat.A,
    #   :A  => mat.A 
      )
end


function calcD(mat::ElasticRod, ips::ElasticRodState)
    return mat.E
end

