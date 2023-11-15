# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


abstract type BC end
@inline Base.:(<<)(a, b::BC) = return (a, b)
@inline Base.:(=>)(a, b::BC) = return (a, b)

# NodeBC


mutable struct NodeBC<:BC
    conds::AbstractDict
    filter::Union{Array{Int,1},Symbol,String,Expr}
    nodes::Array{Node,1}

    function NodeBC(;conds...)
        return new(conds, :(), [])
    end
end


function setup_bc!(model, filter, bc::NodeBC)
    isa(filter, Int) && (filter = [filter])
    bc.filter = filter

    # Filter objects according to bc criteria
    bc.nodes  = model.nodes[bc.filter]
    length(bc.nodes)==0 && notify("setup_bc!: applying boundary conditions to empty array of nodes while evaluating expression ", string(bc.filter))
    # not_found_keys = Set()

    # Find prescribed essential bcs
    for (key,cond) in bc.conds
        for node in bc.nodes
            !haskey(node.dofdict, key) && continue
            dof = node.dofdict[key]
            dof.name==key && (dof.prescribed=true)
        end
    end
end


function compute_bc_vals!(bc::NodeBC, t::Float64, U::Array{Float64,1}, F::Array{Float64,1})
    # essential_keys = Set( dof.name for node in bc.nodes for dof in node.dofs )

    for node in bc.nodes
        x, y, z = node.coord
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


# SurfaceBC and EdgeBC



mutable struct SurfaceBC<:BC
    conds ::AbstractDict
    filter::Union{Symbol,String,Expr}
    faces ::Array{Face,1}

    function SurfaceBC(;conds...)
        return new(conds, :(), [])
    end
end

FaceBC = SurfaceBC


mutable struct EdgeBC<:BC
    conds ::AbstractDict
    filter::Union{Symbol,String,Expr}
    edges ::Array{Edge,1}

    function EdgeBC(;conds...)
        return new(conds, :(), [])
    end
end


function setup_bc!(model, filter, bc::Union{SurfaceBC,EdgeBC})
    bc.filter = filter

    # Filter objects according to bc criteria
    if bc isa SurfaceBC
        bc.faces = model.faces[bc.filter]
        facets = bc.faces
    else
        bc.edges = model.edges[bc.filter]
        facets = bc.edges
    end
    length(facets)==0 && notify("setup_bc!: applying boundary conditions to empty array of faces/edges while evaluating expression ", string(bc.filter))

    # not_found_keys = Set()

    # Find prescribed essential bcs
    for (key,val) in bc.conds
        for facet in facets
            for node in facet.nodes
                !haskey(node.dofdict, key) && continue
                dof = node.dofdict[key]
                dof.name==key && (dof.prescribed=true)
            end
        end
    end
end


function compute_bc_vals!(bc::Union{SurfaceBC,EdgeBC}, t::Float64, U::Array{Float64,1}, F::Array{Float64,1})
    facets = bc isa SurfaceBC ? bc.faces : bc.edges
    essential_keys = Set( dof.name for facet in facets for node in facet.nodes for dof in node.dofs )

    for facet in facets
        elem = facet.owner

        for (key,val) in bc.conds
            if key in essential_keys
                for node in facet.nodes
                    if haskey(node.dofdict, key)
                        dof = node.dofdict[key]
                        x, y, z = node.coord
                        U[dof.eq_id] = eval_arith_expr(val, x=x, y=y, z=z, t=t)
                    end
                end
            else
                Fd, map = distributed_bc(facet, key, val, elem.env, elem.env.ana)
                # Fd, map = distributed_bc(elem, facet, key, val)
                F[map] += Fd
            end
        end
    end
end


# BodyC



mutable struct BodyC<:BC
    conds::AbstractDict
    filter::Union{Array{Int,1},Symbol,String,Expr}
    elems::Array{Element,1}

    function BodyC(;conds...)
        return new(conds, :(), [])
    end
end

ElemBC = BodyC


function setup_bc!(model, filter, bc::BodyC)
    isa(filter, Int) && (filter = [filter])
    bc.filter = filter
    # Filter objects according to bc criteria
    bc.elems = model.elems.active[bc.filter]
    length(bc.elems)==0 && notify("setup_bc!: applying boundary conditions to empty array of elements while evaluating expression ", string(bc.filter))

    # Find prescribed essential bcs
    for (key,val) in bc.conds  # e.g. key = tx, ty, tn, etc...
        for elem in bc.elems
            for node in elem.nodes
                !haskey(node.dofdict, key) && continue # if key not found, assume it is a natural bc
                dof = node.dofdict[key]
                dof.name==key && (dof.prescribed=true)
            end
        end
    end
end


function compute_bc_vals!(bc::BodyC, t::Float64, U::Array{Float64,1}, F::Array{Float64,1})
    essential_keys = Set( dof.name for elem in bc.elems for node in elem.nodes for dof in node.dofs )

    for elem in bc.elems
        for (key,val) in bc.conds
            if key in essential_keys
                for node in elem.nodes
                    if haskey(node.dofdict, key)
                        dof = node.dofdict[key]
                        x, y, z = node.coord
                        U[dof.eq_id] = eval_arith_expr(val, x=x, y=y, z=z, t=t)
                    end
                end
            else
                Fd, map = body_c(elem, key, val)
                F[map] += Fd
            end
        end
    end
end


# Return a vector with all domain dofs and the number of unknown dofs according to bcs
function configure_dofs!(model, bcbinds::Array{<:Tuple,1})

    # get active nodes
    ids = [ node.id for elem in model.elems.active for node in elem.nodes ]
    ids = sort(unique(ids)) # sort is required to preserve node numbering optimization
    active_nodes = model.nodes[ids]

    # All dofs
    # dofs = Dof[dof for node in model.nodes for dof in node.dofs]
    dofs = Dof[dof for node in active_nodes for dof in node.dofs]

    # Reset all dofs as natural conditions
    for dof in dofs
        dof.prescribed = false
    end

    # Setup bcs and prescribed marker for each dof
    for (filter,bc) in bcbinds
        setup_bc!(model, filter, bc)
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
function get_bc_vals(domain, bcs::Array{<:Tuple,1}, t=0.0)
    # This function will be called for each time increment

    ndofs = domain.ndofs
    U = zeros(ndofs)
    F = zeros(ndofs)

    for (key,bc) in bcs
        compute_bc_vals!(bc, t, U, F)
    end

    return U, F
end
