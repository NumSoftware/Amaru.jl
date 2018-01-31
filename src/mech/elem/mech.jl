# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

abstract type Mechanical<:Element end

"""
`elem_config_dofs(elem)`

Sets up the dofs for all nodes in `elem` according with its type.
This function can be specialized by concrete types.
"""
function elem_config_dofs(elem::Mechanical)::Void
    for node in elem.nodes
        add_dof(node, :ux, :fx)
        add_dof(node, :uy, :fy)
        if elem.shared_data.ndim==3; add_dof(node, :uz, :fz) end
    end
end

"""
`elem_init(elem)`

Configures `elem` according to its type.
This function can be specialized by concrete types.
"""
function elem_init(elem::Mechanical)::Void
    # No-op function but can be specialized by concrete types
    # This is called by set_mat(...) function
    return nothing
end

"""
`elem_stiffness(elem)`

Returns the stiffness matrix for `elem`.
This function must be defined by each concrete type.
"""
function elem_stiffness(elem::Mechanical)
    error("elem_stiffness function not defined for material type $(typeof(elem.mat))")
end

"""
`elem_update!(elem, U, F)`

Updates the state of an element given the current global vectors for essential and
natural quantities. It also updates the global vector F.
This function must be defined by each concrete type.
"""
function elem_update!(elem::Mechanical, U::Array{Float64,1}, F::Array{Float64,1}, Δt::Float64)
    error("elem_update function not defined for material type $(typeof(elem.mat))")
end

"""
`elem_RHS(elem)`

Returns the right-hand-side vector for `elem`.
This function can be specialized by concrete types.
"""
function elem_RHS(elem::Mechanical)
    return Float64[], Int64[]
end

"""
`elem_vals(elem)`

Returns a dictionary with values for the element.
Those values are intended to be constant along the element.
This function can be specialized by concrete types.
"""
function elem_vals(elem::Mechanical)
    return Dict{Symbol, Float64}()
end
