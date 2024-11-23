
abstract type AcousticMech<:Element end


"""
`elem_config_dofs(mat, elem)`

Sets up the dofs for all nodes in `elem` according to material mat.
This function can be overloaded by concrete types.
"""
function elem_config_dofs(elem::AcousticMech)
    for node in elem.nodes
        add_dof(node, :up, :fp)
        add_dof(node, :ux, :fx)
        add_dof(node, :uy, :fy)
        elem.ctx.ndim==3 && add_dof(node, :uz, :fz)
    end
end

"""
`elem_init(mat, elem)`

Sets up `elem` according to material `mat`.
This function is called after mat is assigned to `elem` by function `set_mat`.
This function can be overloaded by concrete types.
"""
function elem_init(elem::AcousticMech)
    # No-op function but can be specialized by concrete types
    # This is called by set_mat(...) function
    return nothing
end

"""
`update_elem!(mat, elem)`

Returns the force increment vector dF given a displecement increment vector `dU`
for `elem` according to material `mat`.
This function also updates strains, stresses and internal variables of all
`IpState` objects at integration points.
This function must be redefined by concrete types.
"""
function update_elem!(elem::AcousticMech, dU::Array{Float64,1})
    error("elem_dF function not defined for material type $(typeof(elem.mat))")
end

"""
`elem_vals(mat, elem)`

Returns two dictionaries with keys and values for all nodes and for `elem` itself.
Values for the element are intended to be constant along the element.
This function can be overloaded by concrete types.
"""
function elem_vals(elem::AcousticMech)
    return Dict{Symbol, Float64}()
end


"""
`elem_internal_forces!(elem, F)`

Gets internal nodal forces from current element state.
This function must be defined by each concrete type.
"""
function elem_internal_forces(elem::AcousticMech, F::Array{Float64,1})
end


