# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# Dof
# ===

"""
`Dof()`

Creates an object that represents a Degree of Freedom in a finite element analysis.
`Node` objects include a field called `dofs` which is an array of `Dof` objects.
"""
mutable struct Dof
    name    ::Symbol  # essential variable name
    natname ::Symbol  # natural value name
    eq_id ::Int64     # number of equation in global system
    prescribed::Bool  # flag for prescribed dof
    vals::OrderedDict{Symbol,Float64}
    function Dof(name::Symbol, natname::Symbol)
        new(name, natname, 0, false, OrderedDict{Symbol,Float64}())
    end
end

const NULL_DOF = Dof(:null, :null)

# Node
# ====

"""
`Node(X)`

Creates an object that represents a Node in a finite element analysis. The `X` parameter is a
vector that represents the node coordinates.

**Important fields are**
`X`   : A vector of coordinates
`tag` : An int or string tag
`dofs`: An array of `Dof` objects
"""
mutable struct Node
    id      ::Int
    X       ::Array{Float64,1}
    tag     ::String
    dofs    ::Array{Dof,1}
    dofdict ::OrderedDict{Symbol,Dof}

    function Node(X::Array{Float64,1}; tag::String=0, id::Int=-1)
        this = new(id, X, tag)
        this.dofs = []
        this.dofdict = OrderedDict{Symbol,Dof}()
        return this
    end
    #function Node(point::Point; id::Int=-1)
        #Node([point.x, point.y, point.z], tag=point.tag, id=id)
    #end
end

# The functions below can be used in conjuntion with sort
get_x(node::Node) = node.X[1]
get_y(node::Node) = node.X[2]
get_z(node::Node) = node.X[3]

# Add a new degree of freedom to a node
function add_dof(node::Node, name::Symbol, natname::Symbol)
    if !haskey(node.dofdict, name)
        dof = Dof(name, natname)
        push!(node.dofs, dof)
        node.dofdict[name] = dof
        node.dofdict[natname] = dof
    end
end


# Index operator for node to get a dof
function Base.getindex(node::Node, s::Symbol)
    return node.dofdict[s]
end

# General node sorting
function Base.sort!(nodes::Array{Node,1})
    return sort!(nodes, by=node->sum(node.X))
end


# Get node values in a dictionary
function node_vals(node::Node)
    coords = OrderedDict( :x => node.X[1], :y => node.X[2], :z => node.X[3] )
    all_vals = [ dof.vals for dof in node.dofs ]
    return merge(coords, all_vals...)
end


# Node collection
# ===============

# Index operator for an collection of nodes
function Base.getindex(nodes::Array{Node,1}, s::Symbol)
    s==:all && return nodes
    error("Element getindex: Invalid symbol $s")
end

# Index operator for an collection of nodes
function Base.getindex(nodes::Array{Node,1}, s::String)
    R = [ node for node in nodes if node.tag==s ]
    sort!(R, by=node->sum(node.X))
end

# Index operator for an collection of nodes
function Base.getindex(nodes::Array{Node,1}, filter_ex::Expr)
    R = Node[]
    for node in nodes
        x, y, z = node.X
        eval_arith_expr(filter_ex, x=x, y=y, z=z) && push!(R, node)
    end

    sort!(R, by=node->sum(node.X))
    return R
end

# Get node coordinates for an collection of nodes as a matrix
function nodes_coords(nodes::Array{Node,1}, ndim=3)
    nnodes = length(nodes)
    [ nodes[i].X[j] for i=1:nnodes, j=1:ndim]
end

#=
# Get the dofs ids
@inline function nodes_map(nodes::Array{Node,1}, key::Symbol)
    return [ node.dofdict[key].eq_id for node in elem.nodes if haskey(node.dofdict, key) ]
end

# Get the dofs ids for the given keys
@inline function nodes_map(nodes::Array{Node,1}, keys::NTuple{N, Symbol} ) where N
    return [ node.dofdict[key].eq_id for node in elem.nodes for key in keys if haskey(node.dofdict, key) ]
end

# Get the values for a given key
@inline function nodes_values(nodes::Array{Node,1}, key::Symbol)
    return [ node.dofdict[key].vals[key] for node in elem.nodes if haskey(node.dofdict, s) ]
end

# Get the values for the given keys
@inline function nodes_values(nodes::Array{Node,1}, keys::NTuple{N, Symbol} ) where N
    return [ node.dofdict[key].eq_id for node in elem.nodes for key in keys if haskey(node.dofdict, key) ]
end
=#


function get_data(node::Node)
    table = DTable()
    dict = OrderedDict{Symbol,Float64}(:id=> node.id)
    for dof in node.dofs
        dict = merge(dict, dof.vals)
    end
    push!(table, dict)
    return table
end

function get_data(nodes::Array{Node,1})
    table = DTable()
    for node in nodes
        dict = OrderedDict{Symbol,Float64}(:id=> node.id)
        for dof in node.dofs
            dict = merge(dict, dof.vals)
        end
        push!(table, dict)
    end
    return table
end
