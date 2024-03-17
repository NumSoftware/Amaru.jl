# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# Abstract Element type


"""
    Element

An abstract type to represent a finite element.
Concrete types are specific to certain analysis types.
"""
abstract type Element<:AbstractCell
    #Element subtypes must have the following fields:
    #id    ::Int
    #shape ::CellShape
    #cell  ::Cell
    #nodes ::Array{Node,1}
    #ips   ::Array{Ip,1}
    #tag   ::String
    #mat::Material
    #mat   ::Properties
    #active::Bool
    #linked_elems::Array{Element,1}
    #env ::ModelEnv
end

# @inline Base.:(<<)(a, b::Type{<:Element}) = return (a, b)


# Function to create new concrete types filled with relevant information


function new_element(etype::Type{<:Element}, shape::CellShape, nodes::Array{Node,1}, tag::String, env::ModelEnv)
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


"""
`elem_config_dofs(elem)`

Sets up the dofs for all nodes in `elem` according with its type.
This function can be specialized by concrete types.
"""
function elem_config_dofs(elem::Element)
    return nothing
end


"""
`elem_init(elem)`

Configures `elem` according to its type.
This function can be specialized by concrete types.
"""
function elem_init(elem::Element)
    for ip in elem.ips
        init_state(elem.mat, ip.state)
    end
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


"""
    post_process(elem)

Computes secondary values at the end of each successful increment.
Few special elements need to specialize this function.
"""
function post_process(elem::Element)
    return nothing
end


# Auxiliary functions for elements
# ================================

# Get the element coordinates matrix
function getcoords(elem::Element)
    nnodes = length(elem.nodes)
    ndim   = elem.env.ndim
    return [ elem.nodes[i].coord[j] for i in 1:nnodes, j=1:ndim]
end


function setquadrature!(elem::Element, n::Int=0)

    if !(n in keys(elem.shape.quadrature))
        alert("setquadrature!: cannot set $n integration points for shape $(elem.shape.name)")
        return
    end

    ipc = get_ip_coords(elem.shape, n)
    n = size(ipc,1)

    resize!(elem.ips, n)
    for i in 1:n
        R = ipc[i].coord
        w = ipc[i].w
        elem.ips[i] = Ip(R, w)
        elem.ips[i].id = i
        elem.ips[i].state = compat_state_type(typeof(elem.mat), typeof(elem), elem.env)(elem.env)
        elem.ips[i].owner = elem
    end

    # finding ips global coordinates
    C     = getcoords(elem)
    shape = elem.shape

    # fix for link elements
    if shape.family==LINEJOINTCELL
        bar   = elem.linked_elems[2]
        C     = getcoords(bar)
        shape = bar.shape
    end
    if shape.family==TIPJOINTCELL
        C = reshape(elem.nodes[1].coord, (1,3))
    end

    # fix for joint elements
    if shape.family==JOINTCELL
        C     = C[1:shape.facet_shape.npoints, : ]
        shape = shape.facet_shape
    end

    # interpolation
    for ip in elem.ips
        N = shape.func(ip.R)
        ip.coord = extend!(C'*N, 3)
    end

end


function setquadrature!(elems::Array{<:Element,1}, n::Int=0)
    # shapes = CellShape[]

    for elem in elems
        setquadrature!(elem, n)
        # if n in keys(elem.shape.quadrature)
        # else
        #     if !(elem.shape in shapes)
        #         alert("setquadrature: cannot set $n integration points for shape $(elem.shape.name)")
        #         push!(shapes, elem.shape)
        #     end
        # end
    end
end


function changequadrature!(elems::Array{<:Element,1}, n::Int=0)
    setquadrature!(elems, n)
    foreach(elem_init, elems)
end


function update_material!(elem::Element, mat::Material)
    typeof(elem.mat) == typeof(mat) || error("update_material!: The same material type should be used.")
    elem.mat = mat
end

"""
`changemat(elems, mat)`

Especifies the material model `mat` to be used to represent the behavior of a set of `Element` objects `elems`.
"""
function update_material!(elems::Array{<:Element,1}, mat::Material)
    length(elems)==0 && notify("update_material!: Defining material model $(typeof(mat)) for an empty array of elements.")

    for elem in elems
        update_material!(elem, mat)
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
        fields = fieldnames(compat_state_type(elem.mat))
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
function getnodes(elems::Array{<:Element,1})
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
    s == :ips   && return get_ips(elems)
    return invoke(getproperty, Tuple{Array{<:AbstractCell,1}, Symbol}, elems, s)
end


# General element sorting
function Base.sort!(elems::Array{<:Element,1})
    length(elems)==0 && return

    # General sorting
    sorted = sort(elems, by=elem->sum(getcoords(elem.nodes)))

    # Check type of elements
    shapes = [ elem.shape for elem in elems ]
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

