
abstract type Thermomechanical<:Element end


"""
`elem_config_dofs(mat, elem)`

Sets up the dofs for all nodes in `elem` according to material mat.
This function can be overloaded by concrete types.
"""
function elem_config_dofs(elem::Thermomechanical)
    for node in elem.nodes
        add_dof(node, :ut, :ft)
        add_dof(node, :ux, :fx)
        add_dof(node, :uy, :fy)
        elem.env.ndim==3 && add_dof(node, :uz, :fz)
    end
end

"""
`elem_init(mat, elem)`

Sets up `elem` according to material `mat`.
This function is called after mat is assigned to `elem` by function `set_mat`.
This function can be overloaded by concrete types.
"""
function elem_init(elem::Thermomechanical)
    # No-op function but can be specialized by concrete types
    # This is called by set_mat(...) function
    return nothing
end

"""
`elem_update!(mat, elem)`

Returns the force increment vector dF given a displecement increment vector `dU`
for `elem` according to material `mat`.
This function also updates strains, stresses and internal variables of all
`IpState` objects at integration points.
This function must be redefined by concrete types.
"""
function elem_update!(elem::Thermomechanical, dU::Array{Float64,1})
    error("elem_dF function not defined for material type $(typeof(elem.mat))")
end

"""
`elem_vals(mat, elem)`

Returns two dictionaries with keys and values for all nodes and for `elem` itself.
Values for the element are intended to be constant along the element.
This function can be overloaded by concrete types.
"""
function elem_vals(elem::Thermomechanical)
    return Dict{Symbol, Float64}()
end


function elem_stiffness_matrix(elem::Thermomechanical)
    return zeros(0,0), zeros(Int64,0), zeros(Int64,0)
end


function elem_coupling_matrix(elem::Thermomechanical)
    return zeros(0,0), zeros(Int64,0), zeros(Int64,0)
end


function elem_conductivity_matrix(elem::Thermomechanical)
    return zeros(0,0), zeros(Int64,0), zeros(Int64,0)
end


function elem_RHS_vector(elem::Thermomechanical)
    return zeros(0), zeros(Int64,0)
end
