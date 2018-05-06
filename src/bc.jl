# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

mutable struct BC
    bctype::Symbol
    expr::Union{Symbol,Expr}
    keys::Array{Symbol} 
    vals::Array{T} where T<:Union{Real,Symbol,Expr}
    funs::Array{Functor}

    objects::Array{T} where T<:Union{Node,Face,Edge,Element}

    function BC(bctype::Symbol, target::Union{Int, Array{Int,1}, Expr, Symbol, String}, expconds::Expr=:(); kwconds::Array...)
        bctypes = (:node, :face, :edge, :element)
        !(bctype in bctypes) && error("BC: boundary condition type should be one of $bctypes. :$bctype received")

        ty = typeof(target)
        if ty==Symbol
            targets = (:all, :solids, :lines, :joints, :joints1D)
            !(target in targets) && error("BC: invalid targe :$target")
            expr = target
        elseif ty==Expr
            expr = target
        #elseif ty==Int || ty==Array{Int,1}
            #expr = :(id in $target)
        #elseif ty==String
        elseif ty<:TagType
            expr = :(isequal(tag,$target))
        end

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

        return new(bctype, expr, keys, vals, funs)
    end
end


# Return a vector with all domain dofs and the number of unknown dofs according to bcs
function configure_dofs!(dom, bcs::Array{BC,1})
    # Filter objects according to bc criteria
    for bc in bcs
        dict       = Dict(:node=>dom.nodes, :face=>dom.faces, :edge=>dom.edges, :elem=>dom.elems, :element=>dom.elems)
        bc.objects = dict[bc.bctype][bc.expr]
        length(bc.objects)==0 && warn("eval_bc!: applying boundary conditions to empty array of $(bc.bctype)s while using expression $(bc.expr)")
    end

    # All dofs
    dofs = Dof[dof for node in dom.nodes for dof in node.dofs]

    # Reset all dofs as natural conditions
    for dof in dofs
        dof.prescribed = false
    end

    # Find essential bcs
    for bc in bcs
        if bc.bctype == :node
            for (key,fun) in zip(bc.keys, bc.funs)
                for node in bc.objects
                    !haskey(node.dofdict, key) && error("get_dofs!: key ($key) not found in node $(node.id)")
                    dof = node.dofdict[key]
                    dof.name == key && (dof.prescribed = true)
                end
            end
        else # :face, :edge, :element
            for (key,fun) in zip(bc.keys, bc.funs)  # e.g. key = tx, ty, tn, etc...
                for obj in bc.objects
                    for node in obj.nodes
                        !haskey(node.dofdict, key) && continue # if key not found, assume it is a natural bc
                        dof = node.dofdict[key]
                        dof.name == key && (dof.prescribed = true)
                    end
                end
            end
        end
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
function get_bc_vals(domain, bcs::Array{BC,1}, t=0.0)
    # This function will be called for each time increment

    ndofs = domain.ndofs
    U = zeros(ndofs)
    F = zeros(ndofs)

    essential_keys = Set( dof.name for node in domain.nodes for dof in node.dofs )

    for bc in bcs
        if bc.bctype == :node
            for node in bc.objects
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
        else # bc.type is :face or :edge
            for obj in bc.objects
                
                elem, facet = isa(obj,Element)? (obj,nothing) : (obj.oelem,obj)
                
                for (key,fun) in zip(bc.keys, bc.funs)
                    if key in essential_keys
                        for node in obj.nodes
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
    end

    return U, F
end
