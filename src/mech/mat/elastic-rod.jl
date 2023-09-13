# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export LinearElastic

# """
#     LinearElastic

# A type for linear elastic materials in rods.

# # Fields

# $(TYPEDFIELDS)
# """
# mutable struct LinearElastic<:Material
#     "Young Modulus"
#     E::Float64

#     @doc """
#         $(SIGNATURES)

#     Creates an `LinearElastic` material type

#     # Arguments
#     - `E`: Young modulus
#     - `rho`: Density
#     """
#     function LinearElastic(; params...)
#         names = (E="Young modulus")
#         required = (:E,)
#         @checkmissing params required names

#         params  = values(params)
#         E       = params.E

#         @check E>=0.0
#         return new(E)
#     end
# end






# Type of corresponding state structure
compat_state_type(::Type{LinearElastic}) = ElasticRodState

# Element types that work with this material
compat_elem_types(::Type{LinearElastic}) = (MechRod,)
compat_elem_type_if_embedded(::Type{LinearElastic}) = MechEmbRod


function update_state!(mat::LinearElastic, state::ElasticRodState, Δε::Float64)
    Δσ = mat.E*Δε
    state.ε += Δε
    state.σ += Δσ
    return Δσ, success()
end


function ip_state_vals(mat::LinearElastic, state::ElasticRodState)
    return OrderedDict(
      "sX" => state.σ,
      "eX" => state.ε,
      )
end


function calcD(mat::LinearElastic, ips::ElasticRodState)
    return mat.E
end

