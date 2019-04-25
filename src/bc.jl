# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


abstract type BC end

# NodeBC
# ======

mutable struct NodeBC<:BC
    conds::AbstractDict
    filter::Union{Symbol,String,Expr}
    nodes::Array{Node,1}

    function NodeBC(;conds...)
        return new(conds, :(), [])
    end
end


function setup_bc!(dom, filter, bc::NodeBC)
    bc.filter = filter

    # Filter objects according to bc criteria
    bc.nodes  = dom.nodes[bc.filter]
    length(bc.nodes)==0 && @warn "setup_bc!: applying boundary conditions to empty array of nodes while evaluating expression" bc.filter
    not_found_keys = Set()

    # Find prescribed essential bcs
    for (key,cond) in bc.conds
        for node in bc.nodes
            !haskey(node.dofdict, key) && continue
            dof = node.dofdict[key]
            dof.name == key && (dof.prescribed = true)
        end
    end
end


function compute_bc_vals!(bc::NodeBC, t::Float64, U::Array{Float64,1}, F::Array{Float64,1})
    essential_keys = Set( dof.name for node in bc.nodes for dof in node.dofs )

    for node in bc.nodes
        x, y, z = node.X
        for (key,cond) in bc.conds
            !haskey(node.dofdict, key) && continue
            dof = node.dofdict[key]
            if key==dof.name # essential bc (dof.prescribed should not be modified!)
                U[dof.eq_id] = eval_arith_expr(cond, x=x, y=y, z=z, t=t)
            else # natural bc
                F[dof.eq_id] += eval_arith_expr(cond, x=x, y=y, z=z, t=t)
            end
        end
    end
end


# FaceBC and EdgeBC
# =================


mutable struct FaceBC<:BC
    conds ::AbstractDict
    filter::Union{Symbol,String,Expr}
    faces ::Array{Face,1}

    function FaceBC(;conds...)
        return new(conds, :(), [])
    end
end


mutable struct EdgeBC<:BC
    conds ::AbstractDict
    filter::Union{Symbol,String,Expr}
    edges ::Array{Edge,1}

    function EdgeBC(;conds...)
        return new(conds, :(), [])
    end
end


function setup_bc!(dom, filter, bc::Union{FaceBC,EdgeBC})
    bc.filter = filter

    # Filter objects according to bc criteria
    if bc isa FaceBC
        bc.faces = dom.faces[bc.filter]
        facets = bc.faces
    else
        bc.edges = dom.edges[bc.filter]
        facets = bc.edges
    end
    length(facets)==0 && @warn "setup_bc!: applying boundary conditions to empty array of faces/edges while evaluating expression" bc.filter

    not_found_keys = Set()

    # Find prescribed essential bcs
    for (key,val) in bc.conds
        for facet in facets
            for node in facet.nodes
                !haskey(node.dofdict, key) && continue
                dof = node.dofdict[key]
                dof.name == key && (dof.prescribed = true)
            end
        end
    end
end


function compute_bc_vals!(bc::Union{FaceBC,EdgeBC}, t::Float64, U::Array{Float64,1}, F::Array{Float64,1})
    facets = bc isa FaceBC ? bc.faces : bc.edges
    essential_keys = Set( dof.name for facet in facets for node in facet.nodes for dof in node.dofs )

    for facet in facets
        elem = facet.oelem
        
        for (key,val) in bc.conds
            if key in essential_keys
                for node in facet.nodes
                    if haskey(node.dofdict, key)
                        dof = node.dofdict[key]
                        x, y, z = node.X
                        U[dof.eq_id] = eval_arith_expr(val, x=x, y=y, z=z, t=t)
                    end
                end
            else
                Fd, map = distributed_bc(elem, facet, key, val)
                F[map] += Fd
            end
        end
    end
end


# ElemBC
# ======


mutable struct ElemBC<:BC
    conds::AbstractDict
    filter::Union{Symbol,String,Expr}
    elems::Array{Element,1}

    function ElemBC(;conds...)
        return new(conds, :(), [])
    end
end


# Return a vector with all domain dofs and the number of unknown dofs according to bcs
function setup_bc!(dom, filter, bc::ElemBC)
    bc.filter = filter
    # Filter objects according to bc criteria
    bc.elems = dom.elems[bc.filter]
    length(bc.elems)==0 && @warn "setup_bc!: applying boundary conditions to empty array of elements while evaluating expression" bc.filter

    # Find prescribed essential bcs
    for (key,val) in bc.conds  # e.g. key = tx, ty, tn, etc...
        for elem in bc.elems
            for node in elem.nodes
                !haskey(node.dofdict, key) && continue # if key not found, assume it is a natural bc
                dof = node.dofdict[key]
                dof.name == key && (dof.prescribed = true)
            end
        end
    end
end


function compute_bc_vals!(bc::ElemBC, t::Float64, U::Array{Float64,1}, F::Array{Float64,1})
    essential_keys = Set( dof.name for elem in bc.elems for node in elem.nodes for dof in node.dofs )

    for elem in bc.elems
        for (key,val) in bc.conds
            if key in essential_keys
                for node in elem.nodes
                    if haskey(node.dofdict, key)
                        dof = node.dofdict[key]
                        x, y, z = node.X
                        U[dof.eq_id] = eval_arith_expr(val, x=x, y=y, z=z, t=t)
                    end
                end
            else
                Fd, map = distributed_bc(elem, nothing, key, val)
                F[map] += Fd
            end
        end
    end
end


# Return a vector with all domain dofs and the number of unknown dofs according to bcs
function configure_dofs!(dom, bcbinds::Array{<:Pair,1})

    # All dofs
    dofs = Dof[dof for node in dom.nodes for dof in node.dofs]

    # Reset all dofs as natural conditions
    for dof in dofs
        dof.prescribed = false
    end

    # Setup bcs and prescribed marker for each dof
    for (filter,bc) in bcbinds
        setup_bc!(dom, filter, bc)
    end

    # Split dofs
    presc = [dof.prescribed for dof in dofs]
    pdofs = dofs[presc]
    udofs = dofs[.!presc]
    dofs  = [udofs; pdofs]
    nu    = length(udofs)

    # set eq_id in dofs
    for (i,dof) in enumerate(dofs)
        dof.eq_id = i
    end

    return dofs, nu
end


# Returns the values for essential and natural boundary conditions according to bcs
function get_bc_vals(domain, bcs::Array{<:Pair,1}, t=0.0)
    # This function will be called for each time increment

    ndofs = domain.ndofs
    U = zeros(ndofs)
    F = zeros(ndofs)

    for (key,bc) in bcs
        compute_bc_vals!(bc, t, U, F)
    end

    return U, F
end
