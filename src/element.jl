# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# Abstract Element type
# =====================

abstract type Element<:AbstractCell
    #Element subtypes must have the following fields:
    #id    ::Int
    #shape ::ShapeType
    #cell  ::Cell
    #nodes ::Array{Node,1}
    #ips   ::Array{Ip,1}
    #tag   ::String
    #mat   ::Material
    #active::Bool
    #linked_elems::Array{Element,1}
    #env ::ModelEnv
end

# Function to create new concrete types filled with relevant information

function new_element(etype::Type{<:Element}, shape::ShapeType, nodes::Array{Node,1}, tag::String, env::ModelEnv)
    elem = etype()
    elem.id     = 0
    elem.shape  = shape
    elem.nodes  = nodes
    elem.ips    = []
    elem.tag    = tag
    elem.active = true
    elem.linked_elems = []
    elem.env  = env
    return elem
end


# Functions that should be available in all concrete types derived from Element
# =============================================================================

"""
`elem_config_dofs(elem)`

Sets up the dofs for all nodes in `elem` according with its type.
This function can be specialized by concrete types.
"""
function elem_config_dofs(elem::Element)
    # No-op function but can be specialized by concrete types
    return nothing
end


"""
`elem_init(elem)`

Configures `elem` according to its type.
This function can be specialized by concrete types.
"""
function elem_init(elem::Element)
    # No-op function but can be specialized by concrete types
    return nothing
end


"""
`elem_vals(elem)`

Returns a dictionary with values for the element.
Those values are intended to be constant along the element.
This function can be specialized by concrete types.
"""
function elem_vals(elem::Element)
    return Dict{Symbol, Float64}()
end


"""
`elem_extrapolated_node_vals(elem)`

Returns a dictionary with nodal values obtained by extrapolation
of values at ip points.
"""
function elem_extrapolated_node_vals(elem::Element)
    return Dict{Symbol, Float64}()
end


# Auxiliary functions for elements
# ================================

# Get the element coordinates matrix
function get_coords(elem::Element)
    nnodes = length(elem.nodes)
    ndim   = elem.env.ndim
    return [ elem.nodes[i].coord[j] for i=1:nnodes, j=1:ndim]
end


function setquadrature!(elem::Element, n::Int=0)

    if !(n in keys(elem.shape.quadrature))
        alert("setquadrature!: cannot set $n integration points for shape $(elem.shape.name)")
        return
    end

    ipc = get_ip_coords(elem.shape, n)
    n = size(ipc,1)

    resize!(elem.ips, n)
    for i=1:n
        R = ipc[i,1:3]
        w = ipc[i,4]
        elem.ips[i] = Ip(R, w)
        elem.ips[i].id = i
        elem.ips[i].state = ip_state_type(elem.mat)(elem.env)
        elem.ips[i].owner = elem
    end

    # finding ips global coordinates
    C     = get_coords(elem)
    shape = elem.shape

    # fix for link elements
    if shape.family==JOINT1D_SHAPE
        bar   = elem.linked_elems[2]
        C     = get_coords(bar)
        shape = bar.shape
    end

    # fix for joint elements
    if shape.family==JOINT_SHAPE
        C     = C[1:shape.facet_shape.npoints, : ]
        shape = shape.facet_shape
    end

    # interpolation
    for ip in elem.ips
        N = shape.func(ip.R)
        ip.coord = C'*N
    end

end

function setquadrature!(elems::Array{<:Element,1}, n::Int=0)
    shapes = ShapeType[]

    for elem in elems
        if n in keys(elem.shape.quadrature)
            setquadrature!(elem, n)
        else
            if !(elem.shape in shapes)
                alert("setquadrature: cannot set $n integration points for shape $(elem.shape.name)")
                push!(shapes, elem.shape)
            end
        end
    end
end


function changequadrature!(elems::Array{<:Element,1}, n::Int=0)
    setquadrature!(elems, n)
    foreach(elem_init, elems)
end


function changemat!(elem::Element, mat::Material)
    typeof(elem.mat) == typeof(mat) || error("changemat!: The same material type should be used.")
    elem.mat = mat
end

"""
`changemat(elems, mat)`

Especifies the material model `mat` to be used to represent the behavior of a set of `Element` objects `elems`.
"""
function changemat!(elems::Array{<:Element,1}, mat::Material)
    length(elems)==0 && notify("changemat!: Defining material model $(typeof(mat)) for an empty array of elements.")

    for elem in elems
        changemat!(elem, mat)
    end
end


# Define the state at all integration points in a collection of elements
function setstate!(elems::Array{<:Element,1}; args...)
    length(elems)==0 && notify("setstate!: Setting state to an empty array of elements.")

    greek = Dict(
                 :sigma => :σ,
                 :epsilon => :ε,
                 :eps => :ε,
                 :gamma => :γ,
                 :theta => :θ,
                 :beta => :β,
                )

    found = Set{Symbol}()
    notfound = Set{Symbol}()

    for elem in elems
        fields = fieldnames(ip_state_type(elem.mat))
        for ip in elem.ips
            for (k,v) in args
                if k in fields
                    setfield!(ip.state, k, v)
                    push!(found, k)
                else
                    gk = get(greek, k, :none)
                    if gk == :none
                        push!(notfound, k)
                    else
                        if gk in fields
                            setfield!(ip.state, gk, v)
                            push!(found, k)
                        else
                            push!(notfound, k)
                        end
                    end
                end
            end
        end
    end

    for k in notfound
        if k in found
            msg1 = k in keys(greek) ? " ($(greek[k])) " : ""
            alert("setstate!: field '$k$msg1' was not found at some elements while setting state values")
        else
            error("setstate!: field '$k$msg1' was not found at selected elements while setting state values")
        end
    end
end


# Get all nodes from a collection of elements
function get_nodes(elems::Array{<:Element,1})
    nodes = Set{Node}()
    for elem in elems
        for node in elem.nodes
            push!(nodes, node)
        end
    end
    return [node for node in nodes]
end

# Get all dofs from an element
function get_dofs(elem::Element)
    return Dof[ dof for node in elem.nodes for dof in node.dofs ]
end

function get_dofs(elem::Element, dofname::Symbol)
    return Dof[ dof for node in elem.nodes for dof in node.dofs if dof.name==dofname ]
end


# Get all ips from a collection of elements
function get_ips(elems::Array{<:Element,1})
    return Ip[ ip for elem in elems for ip in elem.ips ]
end


function Base.getproperty(elems::Array{<:Element,1}, s::Symbol)
    s == :solids && return filter(elem -> elem.shape.family==SOLID_SHAPE, elems)
    s == :lines && return filter(elem -> elem.shape.family==LINE_SHAPE, elems)
    s == :embedded && return filter(elem -> elem.shape.family==LINE_SHAPE && length(elem.linked_elems)>0, elems)
    s in (:joints1d, :joints1D) && return filter(elem -> elem.shape.family==JOINT1D_SHAPE, elems)
    s == :joints && return filter(elem -> elem.shape.family==JOINT_SHAPE, elems)
    s == :nodes && return get_nodes(elems)
    s == :ips   && return get_ips(elems)
    error("type Array{Element,1} has no property $s")
end


# Index operator for a collection of elements
function Base.getindex(elems::Array{<:Element,1}, s::Symbol)
    s == :all && return elems
    s == :solids && return filter(elem -> elem.shape.family==SOLID_SHAPE, elems)
    s == :lines && return filter(elem -> elem.shape.family==LINE_SHAPE, elems)
    s == :embedded && return filter(elem -> elem.shape.family==LINE_SHAPE && length(elem.linked_elems)>0, elems)
    s in (:joints1d, :joints1D) && return filter(elem -> elem.shape.family==JOINT1D_SHAPE, elems)
    s == :joints && return filter(elem -> elem.shape.family==JOINT_SHAPE, elems)
    s == :nodes && return get_nodes(elems)
    s == :ips && return get_ips(elems)
    error("Element getindex: Invalid symbol $s")
end


# Index operator for a collection of elements using a string
function Base.getindex(elems::Array{<:Element,1}, tag::String)
    return [ elem for elem in elems if elem.tag==tag ]
end


# Index operator for a collection of elements using an expression
function Base.getindex(elems::Array{<:Element,1}, filter_ex::Expr)
    length(elems)==0 && return Element[]

    nodes = elems[:nodes]
    nodemap = zeros(Int, maximum(node.id for node in nodes) )
    T = Bool[]
    for (i,node) in enumerate(nodes)
        nodemap[node.id] = i
        x, y, z = node.coord
        push!(T, eval_arith_expr(filter_ex, x=x, y=y, z=z))
    end

    R = Element[]
    for elem in elems
        all( T[nodemap[node.id]] for node in elem.nodes ) && push!(R, elem)
    end
    return R
end

# General element sorting
function Base.sort!(elems::Array{<:Element,1})
    length(elems)==0 && return

    # General sorting
    sorted = sort(elems, by=elem->sum(get_coords(elem.nodes)))

    # Check type of elements
    shapes = [ elem.shapetype for elem in elems ]
    shape  = shapes[1]

    if all(shapes.==shape)
        if shape in (LIN2, LIN3)
            node_ids = Set(node.id for elem in sorted for node in elem.nodes)
            for elem in sorted

            end
        elseif shape in (JLIN3, JLIN4)
        end
    end

    # General sorting
    return sorted
end


"""
`elem_ip_vals(elem)`

Returns a table with values from all ips in the element.
This function can be specialized by concrete types.
"""
function elems_ip_vals(elem::Element)
    table = DataTable()
    for ip in elem.ips
        D = ip_state_vals(elem.mat, ip.state)
        push!(table, D)
    end

    return table
end


