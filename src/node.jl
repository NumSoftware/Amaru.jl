# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# Dof
# ===

"""
`Dof()`

Creates an object that represents a Degree of Freedom in a finite element analysis.
`Node` objects include a field called `dofs` which is an array of `Dof` objects.
"""
mutable struct Dof
    name      ::Symbol # essential variable name
    natname   ::Symbol # natural value name
    eq_id     ::Int64  # number of equation in global system
    prescribed::Bool   # flag for prescribed dof
    vals      ::OrderedDict{Symbol,Float64}
    function Dof(name::Symbol, natname::Symbol)
        new(name, natname, 0, false, OrderedDict{Symbol,Float64}())
    end
end


function Base.copy(dof::Dof)
    newdof = Dof(dof.name, dof.natname)
    newdof.eq_id = dof.eq_id
    newdof.prescribed = dof.prescribed
    newdof.vals = copy(dof.vals)
    return newdof
end


function Base.getindex(dofs::Array{Dof,1}, s::Symbol)
    for dof in dofs
        dof.name == s && return dof
        dof.natname == s && return dof
    end
    error("getindex: Dof key $s not found.")
end


function Base.haskey(dofs::Array{Dof,1}, s::Symbol)
    for dof in dofs
        dof.name == s && return true
        dof.natname == s && return true
    end
    return false
end


const null_Dof = Dof(:null, :null)
@inline null(::Type{Dof}) = null_Dof


"""
    Node

A type that represents a finite element node.

# Fields
$(TYPEDFIELDS)
"""
mutable struct Node<:AbstractPoint
    "identification number"
    id      ::Int
    "coordinates vector"
    coord   ::Vec3
    "string tag used to group nodes"
    tag     ::String
    "array of degrees of freedom"
    dofs    ::Array{Dof,1}
    "dictionary of degrees of freedom"
    dofdict ::OrderedDict{Symbol,Dof}

    @doc """
        $(TYPEDSIGNATURES)

    Constructs an uninitiallised `Node`.
    """
    function Node()
        this = new()
        this.id = -1
        this.coord = Vec3()
        this.dofs = Dof[]
        this.dofdict = OrderedDict{Symbol,Dof}()
        return this
    end

    @doc """
        $(TYPEDSIGNATURES)

    Constructs an `Node` with coordinates `x`, `y` and `z`. 
    A `tag` string can be provided optionally.

    # Examples
    
    ```jldoctest
    julia> using Amaru;
    julia> Node(1.0, 2.0, 3.0, tag="vertex")
    Node
      id: -1
      coord: [1.0, 2.0, 3.0]
      tag: "vertex"
      dofs: 0-element Vector{Dof}
      dofdict: OrderedDict{Symbol, Dof} with 0 entries
    ```
    """
    function Node(x::Real, y::Real=0.0, z::Real=0.0; tag::String="", id::Int=-1)
        x = round(x, digits=8) + 0.0 # +0.0 required to drop negative bit
        y = round(y, digits=8) + 0.0
        z = round(z, digits=8) + 0.0
        
        this = new(id, Vec3(x,y,z), tag)
        this.dofs = Dof[]
        this.dofdict = OrderedDict{Symbol,Dof}()
        return this
    end

    @doc """
        $(TYPEDSIGNATURES)

    Constructs an `Node` with coordinates privided in the `X` vector. 
    A `tag` string can be provided optionally.

    # Examples
    
    ```jldoctest
    julia> using Amaru;
    julia> Node([1.0, 2.0, 3.0], tag="vertex")
    Node
      id: -1
      coord: [1.0, 2.0, 3.0]
      tag: "vertex"
      dofs: 0-element Vector{Dof}
      dofdict: OrderedDict{Symbol, Dof} with 0 entries
    ```
    """
    function Node(X::AbstractArray; tag::String="", id::Int=-1)
        @assert length(X) in (1,2,3)
        return Node(X...; tag=tag, id=id)
    end
end


const null_Node = Node(NaN, NaN, NaN)
@inline null(::Type{Node}) = null_Node


#Base.hash(n::Node) = hash( (round(n.coord.x, digits=8), round(n.coord.y, digits=8), round(n.coord.z, digits=8)) )
# Base.hash(n::Node) = hash( (n.coord.x, n.coord.y, n.coord.z) )
Base.hash(n::Node) = hash( (n.coord.x+1.0, n.coord.y+2.0, n.coord.z+3.0) ) # 1,2,3 aim to avoid clash in some arrays of nodes.

"""
    $(TYPEDSIGNATURES)

Creates a copy of `node`.
"""
function Base.copy(node::Node)
    newnode = Node(node.coord, tag=node.tag, id=node.id)
    for dof in node.dofs
        newdof = copy(dof)
        push!(newnode.dofs, newdof)
        newnode.dofdict[dof.name] = newdof
        newnode.dofdict[dof.natname] = newdof
    end
    return newnode
end


# The functions below can be used in conjuntion with sort
get_x(node::Node) = node.coord[1]
get_y(node::Node) = node.coord[2]
get_z(node::Node) = node.coord[3]

"""
    $(TYPEDSIGNATURES)

Adds to `node` a new degree of freedom called `name`
with associated natural variable `natname`.
"""
function add_dof(node::Node, name::Symbol, natname::Symbol)
    if !haskey(node.dofdict, name)
        dof = Dof(name, natname)
        push!(node.dofs, dof)
        node.dofdict[name] = dof
        node.dofdict[natname] = dof
    end
    return nothing
end


# Index operator for node to get a dof
function Base.getindex(node::Node, s::Symbol)
    return node.dofdict[s]
end

# General node sorting
function Base.sort!(nodes::Array{Node,1})
    return sort!(nodes, by=node->sum(node.coord))
end


# Get node values in a dictionary
function node_vals(node::Node)
    coords = OrderedDict( :x => node.coord[1], :y => node.coord[2], :z => node.coord[3] )
    all_vals = [ dof.vals for dof in node.dofs ]
    return merge(coords, all_vals...)
end


# Node collection
# ===============

# Index operator for an collection of nodes
#function Base.getindex(nodes::Array{Node,1}, s::Symbol)
    #s==:all && return nodes
    #error("Element getindex: Invalid symbol $s")
#end

# Index operator for an collection of nodes
function Base.getindex(nodes::Array{Node,1}, s::String)
    R = [ node for node in nodes if node.tag==s ]
    sort!(R, by=node->sum(node.coord))
end

# Index operator for an collection of nodes
function Base.getindex(nodes::Array{Node,1}, filter_ex::Expr)
    R = Node[]
    for node in nodes
        x, y, z = node.coord
        eval_arith_expr(filter_ex, x=x, y=y, z=z) && push!(R, node)
    end

    sort!(R, by=node->sum(node.coord))
    return R
end

# Get node coordinates for a collection of nodes as a matrix
function getcoords(nodes::Array{Node,1}, ndim=3)
    nnodes = length(nodes)
    [ nodes[i].coord[j] for i=1:nnodes, j=1:ndim]
end

function setcoords!(nodes::Array{Node,1}, coords::AbstractArray{Float64,2})
    nrows, ncols = size(coords)
    @assert nrows == length(nodes)

    for (i,node) in enumerate(nodes)
        node.coord.x = coords[i,1]
        node.coord.y = coords[i,2]
        ncols==3 && (node.coord.z = coords[i,3])
    end
end

function get_data(node::Node)
    table = DataTable()
    dict = OrderedDict{Symbol,Real}(:id => node.id)
    for dof in node.dofs
        dict = merge(dict, dof.vals)
    end
    push!(table, dict)
    return table
end

function get_data(nodes::Array{Node,1})
    table = DataTable()
    for node in nodes
        dict = OrderedDict{Symbol,Real}(:id => node.id)
        for dof in node.dofs
            dict = merge(dict, dof.vals)
        end
        push!(table, dict)
    end
    return table
end


function Base.getproperty(nodes::Array{Node,1}, s::Symbol)
    s == :dofs && return [ dof for node in nodes for dof in node.dofs ]
    error("type Array{Node,1} has no property $s")
end

function setvalue!(dof::Dof, sym_val::Pair)
    sym, val = sym_val
    if haskey(dof.vals, sym)
        dof.vals[sym] = val
    end
end

function setvalue!(dofs::Array{Dof,1}, sym_val::Pair)
    for dof in dofs
        setvalue!(dof, sym_val)
    end
end