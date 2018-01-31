# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

import Base.reset
import Base.getindex
import Base.maximum
import Base.minimum
import Base.sort

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
    vals::Dict{Symbol,Float64}
    function Dof(name::Symbol, natname::Symbol) 
        new(name, natname, 0, false, Dict())
    end
end


# Node
# ====

"""
`Node(X)` 

Creates an object that represents a Node in a finite element analysis. The `X` parameter is a 
vector that represents the node coordinates.

**Important fields are**
`X`   : A vector of coordinates
`tag` : A string tag
`dofs`: An array of `Dof` objects
"""
mutable struct Node
    X       ::Array{Float64,1}
    tag     ::AbstractString
    id      ::Int
    dofs    ::Array{Dof,1}
    dofdict ::Dict{Symbol,Dof}
 
    function Node(X::Array{Float64,1}; tag::String="", id::Int=-1)
        this = new(X, tag, id)
        this.dofs = []
        this.dofdict = Dict()
        return this
    end
    #function Node(point::Point; id::Int=-1)
        #Node([point.x, point.y, point.z], tag=point.tag, id=id)
    #end
end

# The functions below can be used in conjuntion with sort
import FemMesh.get_x, FemMesh.get_y, FemMesh.get_z
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
function getindex(node::Node, s::Symbol)
    return node.dofdict[s]
end


# Get node values in a dictionary
function node_vals(node::Node)
    coords = Dict( :x => node.X[1], :y => node.X[2], :z => node.X[3] )
    all_vals = [ dof.vals for dof in node.dofs ]
    return merge(coords, all_vals...)
end


# Node collection
# ===============

# Index operator for an collection of nodes
function getindex(nodes::Array{Node,1}, s::Symbol) 
    s==:all && return nodes
    error("Element getindex: Invalid symbol $s")
end

# Index operator for an collection of nodes
function getindex(nodes::Array{Node,1}, filter_ex::Expr) 
    @assert filter_ex.head in (:call, :&&, :||)
    expr = fix_comparison_scalar(filter_ex)
    fun  = Functor(:(x,y,z,id,tag), expr)
    return [ n for n in nodes if fun(n.X[1], n.X[2], n.X[3], n.id, n.tag) ]
end

# Get node coordinates for an collection of nodes as a matrix
function nodes_coords(nodes::Array{Node,1}, ndim=3)
    nnodes = length(nodes)
    [ nodes[i].X[j] for i=1:nnodes, j=1:ndim]
end
