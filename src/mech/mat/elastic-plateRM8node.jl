# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# Rectangular Reissner Mindlin Plate FEM

export ElasticPlateRM8node

mutable struct ElasticPlateRM8nodeState<:IpState
    env::ModelEnv
    function ElasticPlateRM8nodeState(env::ModelEnv=ModelEnv())
        return new(env)
    end
end

mutable struct ElasticPlateRM8node<:Material
    E::Float64
    nu::Float64
    ρ::Float64

    function ElasticPlateRM8node(prms::Dict{Symbol,Float64})
        return  ElasticPlateRM8node(;prms...)
    end

    function ElasticPlateRM8node(;E=NaN, nu=NaN, ρ=0.0)
        E>0.0 || error("Invalid value for E: $E")
        (0<=nu<0.5) || error("Invalid value for nu: $nu")

        this = new(E, nu, ρ)
        return this
    end
end

matching_elem_type(::ElasticPlateRM8node, shape::CellShape, ndim::Int) = PlateRM8node

# Type of corresponding state structure
ip_state_type(mat::ElasticPlateRM8node) = ElasticPlateRM8nodeState

function ip_state_vals(mat::ElasticPlateRM8node, state::ElasticPlateRM8nodeState)
    return OrderedDict{Symbol, Float64}()
end
