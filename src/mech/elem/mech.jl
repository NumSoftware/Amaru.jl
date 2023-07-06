# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

abstract type MechElem<:Element end


"""
`elem_config_dofs(elem)`

Sets up the dofs for all nodes in `elem` according with its type.
This function can be specialized by concrete types.
"""
function elem_config_dofs(elem::MechElem)
    for node in elem.nodes
        add_dof(node, :ux, :fx)
        elem.env.ndim>=2 && add_dof(node, :uy, :fy)
        elem.env.ndim==3 && add_dof(node, :uz, :fz)
    end
end


"""
`elem_stiffness(elem)`

Returns the stiffness matrix for `elem`.
This function must be defined by each concrete type.
"""
function elem_stiffness(elem::MechElem)
    error("elem_stiffness function not defined for material type $(typeof(elem.matparams))")
end

"""
`elem_mass(elem)`

Returns the mass matrix for `elem`.
This function must be defined by each concrete type.
"""
function elem_mass(elem::MechElem)
   ndim=elem.env.ndim
   ndofs = length(elem.nodes)*ndim
   M = zeros(ndofs, ndofs)
   keys = (:ux, :uy, :uz)[1:ndim]
   map  = [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
   return M, map, map
end

"""
`elem_internal_forces!(elem, F)`

Gets internal nodal forces from current element state.
This function must be defined by each concrete type.
"""
function elem_internal_forces(elem::MechElem, F::Array{Float64,1})
end

"""
`update_elem!(elem, U, F)`

Updates the state of an element given the current global vectors for essential and
natural quantities. It also updates the global vector F.
This function must be defined by each concrete type.
"""
function update_elem!(elem::MechElem, U::Array{Float64,1}, Δt::Float64)
    error("update_elem function not defined for material type $(typeof(elem.matparams))")
end
