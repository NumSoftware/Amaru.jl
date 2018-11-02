# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru


function parse_conditions(expconds::Expr=:(); kwconds...)
    # List of boundary conditions
    keys = Symbol[]
    vals = Union{Real,Symbol,Expr}[]
    funs = Functor[]

    conds = expconds.head == :(=) ? Expr(:tuple, expconds) : expconds
    conds.head != :tuple && error("BC: invalid expression in boundary condition: $expconds")

    # evaluates each condition in tuple
    for cond in conds.args
        cond.head != :(=) && error("BC: invalid expression in boundary condition: $cond")
        key, val = cond.args[1:2] 
        typeof(key) != Symbol && error("BC: invalid expression in boundary condition: $cond")
        fun = Functor( :(t,x,y,z), val )
        push!(keys, key)
        push!(vals, val)
        push!(funs, fun)
    end

    # add conditions from keyword arguments
    for (key,val) in kwconds
        fun = Functor( :(t,x,y,z), val )
        push!(keys, key)
        push!(vals, val)
        push!(funs, fun)
    end

    return keys, vals, funs
end


abstract type BC
end


# NodeBC
# ======


mutable struct NodeBC<:BC
    expr::Union{Symbol,Expr}
    keys::Array{Symbol} 
    vals::Array{T} where T<:Union{Real,Symbol,Expr}
    funs::Array{Functor}
    nodes::Array{Node,1}

    function NodeBC(expr::Union{Expr,TagType}, expconds::Expr=:(); kwconds...)
        typeof(expr)<:TagType && (expr=:(isequal(tag,$expr)))
        keys, vals, funs = parse_conditions(expconds; kwconds...)
        return new(expr, keys, vals, funs, [])
    end

    function NodeBC(nodes::Array{Node,1}, expconds::Expr=:(); kwconds...)
        keys, vals, funs = parse_conditions(expconds; kwconds...)
        return new(:(), keys, vals, funs, nodes)
    end
end


function setup_bc!(dom, bc::NodeBC)
    # Filter objects according to bc criteria
    bc.expr==:() || (bc.nodes = dom.nodes[bc.expr])
    length(bc.nodes)==0 && @warn "setup_bc!: applying boundary conditions to empty array of nodes while evaluating expression" bc.expr

    # Find prescribed essential bcs
    for (key,fun) in zip(bc.keys, bc.funs)
        for node in bc.nodes
            !haskey(node.dofdict, key) && continue #error("get_dofs!: key ($key) not found in node $(node.id)")
            dof = node.dofdict[key]
            dof.name == key && (dof.prescribed = true)
        end
    end
end


function compute_bc_vals!(bc::NodeBC, t::Float64, U::Array{Float64,1}, F::Array{Float64,1})
    essential_keys = Set( dof.name for node in bc.nodes for dof in node.dofs )

    for node in bc.nodes
        x, y, z = node.X
        for (key,fun) in zip(bc.keys, bc.funs)
            dof = node.dofdict[key]
            if key==dof.name # essential bc (dof.prescribed should not be modified!)
                U[dof.eq_id] = fun(t,x,y,z)
            else # natural bc
                F[dof.eq_id] += fun(t,x,y,z)
            end
        end
    end
end


# FaceBC and EdgeBC
# =================


mutable struct FaceBC<:BC
    expr::Union{Symbol,Expr}
    keys::Array{Symbol} 
    vals::Array{T} where T<:Union{Real,Symbol,Expr}
    funs::Array{Functor}
    faces::Array{Face,1}

    function FaceBC(expr::Union{Expr,TagType}, expconds::Expr=:(); kwconds...)
        typeof(expr)<:TagType && (expr=:(isequal(tag,$expr)))
        keys, vals, funs = parse_conditions(expconds; kwconds...)
        return new(expr, keys, vals, funs, [])
    end

    function FaceBC(faces::Array{<:Face,1}, expconds::Expr=:(); kwconds...)
        keys, vals, funs = parse_conditions(expconds; kwconds...)
        return new(:(), keys, vals, funs, faces)
    end
end


mutable struct EdgeBC<:BC
    expr::Union{Symbol,Expr}
    keys::Array{Symbol} 
    vals::Array{T} where T<:Union{Real,Symbol,Expr}
    funs::Array{Functor}
    edges::Array{Edge,1}

    function EdgeBC(expr::Union{Expr,TagType}, expconds::Expr=:(); kwconds...)
        typeof(expr)<:TagType && (expr=:(isequal(tag,$expr)))
        keys, vals, funs = parse_conditions(expconds; kwconds...)
        return new(expr, keys, vals, funs, [])
    end

    function EdgeBC(edges::Array{<:Edge,1}, expconds::Expr=:(); kwconds...)
        keys, vals, funs = parse_conditions(expconds; kwconds...)
        return new(:(), keys, vals, funs, edges)
    end
end


function setup_bc!(dom, bc::Union{FaceBC,EdgeBC})
    # Filter objects according to bc criteria
    if typeof(bc)==FaceBC
        bc.expr==:() || (bc.faces = dom.faces[bc.expr])
        facets = bc.faces
    else
        bc.expr==:() || (bc.edges = dom.edges[bc.expr])
        facets = bc.edges
    end
    length(facets)==0 && @warn "setup_bc!: applying boundary conditions to empty array of faces/edges while evaluating expression" bc.expr

    #@show bc.expr
    #@show facets

    # Find prescribed essential bcs
    for (key,fun) in zip(bc.keys, bc.funs)  # e.g. key = tx, ty, tn, etc...
        for facet in facets
            for node in facet.nodes
                !haskey(node.dofdict, key) && continue # if key not found, assume it is a natural bc
                dof = node.dofdict[key]
                dof.name == key && (dof.prescribed = true)
            end
        end
    end
end


function compute_bc_vals!(bc::Union{FaceBC,EdgeBC}, t::Float64, U::Array{Float64,1}, F::Array{Float64,1})
    facets = typeof(bc)==FaceBC ? bc.faces : bc.edges
    essential_keys = Set( dof.name for facet in facets for node in facet.nodes for dof in node.dofs )

    for facet in facets
        elem = facet.oelem
        
        for (key,fun) in zip(bc.keys, bc.funs)
            if key in essential_keys
                for node in facet.nodes
                    if haskey(node.dofdict, key)
                        dof = node.dofdict[key]
                        x, y, z = node.X
                        U[dof.eq_id] = fun(t,x,y,z)
                    end
                end
            else
                Fd, map = distributed_bc(elem, facet, key, fun)
                F[map] += Fd
            end
        end
    end
end


# ElemBC
# ======


mutable struct ElemBC<:BC
    expr::Union{Symbol,Expr}
    keys::Array{Symbol} 
    vals::Array{T} where T<:Union{Real,Symbol,Expr}
    funs::Array{Functor}
    elems::Array{Element,1}

    function ElemBC(expr::Union{Expr,TagType}, expconds::Expr=:(); kwconds...)
        typeof(expr)<:TagType && (expr=:(isequal(tag,$expr)))
        keys, vals, funs = parse_conditions(expconds; kwconds...)
        return new(expr, keys, vals, funs, [])
    end

    function ElemBC(elems::Array{<:Element,1}, expconds::Expr=:(); kwconds...)
        keys, vals, funs = parse_conditions(expconds; kwconds...)
        return new(:(), keys, vals, funs, elems)
    end
end


# Return a vector with all domain dofs and the number of unknown dofs according to bcs
function setup_bc!(dom, bc::ElemBC)
    # Filter objects according to bc criteria
    bc.expr==:() || (bc.elems = dom.elems[bc.expr])
    length(bc.elems)==0 && @warn "setup_bc!: applying boundary conditions to empty array of elements while evaluating expression" bc.expr

    # Find prescribed essential bcs
    for (key,fun) in zip(bc.keys, bc.funs)  # e.g. key = tx, ty, tn, etc...
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
        for (key,fun) in zip(bc.keys, bc.funs)
            if key in essential_keys
                for node in elem.nodes
                    if haskey(node.dofdict, key)
                        dof = node.dofdict[key]
                        x, y, z = node.X
                        U[dof.eq_id] = fun(t,x,y,z)
                    end
                end
            else
                Fd, map = distributed_bc(elem, nothing, key, fun)
                F[map] += Fd
            end
        end
    end
end


# Return a vector with all domain dofs and the number of unknown dofs according to bcs
function configure_dofs!(dom, bcs::Array{<:BC,1})

    # All dofs
    dofs = Dof[dof for node in dom.nodes for dof in node.dofs]

    # Reset all dofs as natural conditions
    for dof in dofs
        dof.prescribed = false
    end

    # Setup bcs and prescribed marker for each dof
    for bc in bcs
        setup_bc!(dom, bc)
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
function get_bc_vals(domain, bcs::Array{<:BC,1}, t=0.0)
    # This function will be called for each time increment

    ndofs = domain.ndofs
    U = zeros(ndofs)
    F = zeros(ndofs)

    for bc in bcs
        compute_bc_vals!(bc, t, U, F)
    end

    return U, F
end
