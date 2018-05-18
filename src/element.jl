# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# Abstract Element type
# =====================

abstract type Element 
    #Element subtypes must have the following fields:
    #id    ::Int
    #shape ::ShapeType
    #nodes ::Array{Node,1}
    #ips   ::Array{Ip,1}
    #tag   ::String
    #mat   ::Material
    #active::Bool
    #linked_elems::Array{Element,1}
    #shared_data ::SharedAnalysisData
end

# Function to create new concrete types filled with relevant information
function new_element(etype::Type{<:Element}, shape::ShapeType, nodes::Array{Node,1}, shared_data::SharedAnalysisData, tag::TagType=0)
    elem = etype()
    elem.id     = 0
    elem.shape  = shape
    elem.nodes  = nodes
    elem.ips    = []
    elem.tag    = tag
    elem.active = true
    elem.linked_elems = []
    elem.shared_data  = shared_data
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
function elem_coords(elem::Element)
    nnodes = length(elem.nodes)
    ndim   = elem.shared_data.ndim
    return [ elem.nodes[i].X[j] for i=1:nnodes, j=1:ndim]
end


function elem_config_ips(elem::Element, nips::Int=0)
    ipc =  get_ip_coords(elem.shape, nips)
    nips = size(ipc,1)

    resize!(elem.ips, nips)
    for i=1:nips
        R = ipc[i,1:3]
        w = ipc[i,4]
        elem.ips[i] = Ip(R, w)
        elem.ips[i].id = i
        elem.ips[i].data = new_ip_state(elem.mat, elem.shared_data)
        elem.ips[i].owner = elem
    end

    # finding ips global coordinates
    C     = elem_coords(elem)
    shape = elem.shape

    # fix for link elements
    if shape.class==JOINT1D_SHAPE
        bar   = elem.linked_elems[2]
        C     = elem_coords(bar)
        shape = bar.shape
    end

    # fix for joint elements
    if shape.class==JOINT_SHAPE
        C     = C[1:div(end,2),:]
        shape = shape.facet_shape
    end 

    # interpolation
    for ip in elem.ips
        N = shape.func(ip.R)
        ip.X = C'*N
        if length(ip.X)==2 # complete z=0.0 for 2D analyses
            ip.X = [ ip.X; 0.0 ]
        end
    end
end


"""
`set_mat(elems, mat, [nips=0])`

Especifies the material model `mat` to be used to represent the behavior of a set of `Element` objects `elems`.
"""
function set_mat(elems::Array{Element,1}, mat::Material; nips::Int64=0)
    if length(elems)==0
        warn("Defining material model $(typeof(mat)) for an empty array of elements.\n")
    end
    for elem in elems
        set_mat(elem, mat, nips=nips)
    end
end


# Define the state at all integration points in a collection of elements
function set_state(elems::Array{Element,1}; args...)
    for elem in elems
        for ip in elem.ips
            set_state(ip.data; args...)
        end
    end
end


# Get all nodes from a collection of elements
function get_nodes(elems::Array{Element,1})
    nodes = Set{Node}()
    for elem in elems
        for node in elem.nodes
            push!(nodes, node)
        end
    end
    return [node for node in nodes]
end


# Get all ips from a collection of elements
function get_ips(elems::Array{Element,1})
    ips = Ip[]
    for elem in elems
        for ip in elem.ips
            push!(ips, ip)
        end
    end
    return ips
end

# Index operator for an element
function getindex(elem::Element, s::Symbol)
    if s == :nodes
        return elem.nodes
    end
    if s == :ips
        return elem.ips
    end
    error("Element getindex: Invalid symbol $s")
end

# Index operator for a collection of elements
function getindex(elems::Array{Element,1}, s::Symbol)
    s == :all && return elems
    s == :solids && return filter(elem -> elem.shape.class==SOLID_SHAPE, elems)
    s == :lines && return filter(elem -> elem.shape.class==LINE_SHAPE, elems)
    s == :embedded && return filter(elem -> elem.shape.class==LINE_SHAPE && length(elem.linked_elems)>0, elems)
    s in (:joints1d, :joints1D) && return filter(elem -> elem.shape.class==JOINT1D_SHAPE, elems)
    s == :joints && return filter(elem -> elem.shape.class==JOINT_SHAPE, elems)
    s == :nodes && return get_nodes(elems)
    s == :ips && return get_ips(elems)
    error("Element getindex: Invalid symbol $s")
end

# Index operator for a collection of elements using an expression
function getindex(elems::Array{Element,1}, filter_ex::Expr)
    @assert filter_ex.head in (:call, :&&, :||)
    expr = fix_comparison_arrays(filter_ex)
    fun  = Functor(:(x,y,z,id,tag), expr)

    result = Array{Element}(0)
    for elem in elems
        coords = nodes_coords(elem.nodes)
        x = coords[:,1]
        y = coords[:,2]
        z = coords[:,3]
        fun(x, y, z, elem.id, elem.tag) && push!(result, elem) 
    end
    return result
end

"""
`elem_ip_vals(elem)`

Returns a table with values from all ips in the element.
This function can be specialized by concrete types.
"""
function elems_ip_vals(elem::Element)
    table = DTable()
    for ip in elem.ips
        D = ip_state_vals(elem.mat, ip.data)
        push!(table, D)
    end

    return table
end

