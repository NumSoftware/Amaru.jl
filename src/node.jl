# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

"""
    Node

A type that represents a finite element node.

# Fields
$(TYPEDFIELDS)
"""
mutable struct Node<:AbstractPoint
    id     ::Int
    coord  ::Vec3
    tag    ::String
    dofs   ::Array{Dof,1}
    dofdict::OrderedDict{Symbol,Dof}
    vals   ::OrderedDict{Symbol,Float64}
    aux    ::Bool
    # elems  ::Array{AbstractCell,1} # "elements that share the node"

    @doc """
        $(TYPEDSIGNATURES)

    Constructs an uninitiallised `Node`.
    """
    function Node()
        this         = new()
        this.id      = -1
        this.coord   = Vec3()
        this.dofs    = Dof[]
        this.dofdict = OrderedDict{Symbol,Dof}()
        this.vals    = OrderedDict{Symbol,Float64}()
        this.aux     = false
        # this.elems   = []
        return this
    end

    @doc """
        $(TYPEDSIGNATURES)

    Constructs a `Node` with coordinates `x`, `y` and `z`.
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
        
        this         = new(id, Vec3(x,y,z), tag)
        this.dofs    = Dof[]
        this.dofdict = OrderedDict{Symbol,Dof}()
        this.aux     = false
        return this
    end

    @doc """
        $(TYPEDSIGNATURES)

    Constructs a `Node` with coordinates provided in the `X` vector.
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
Base.isequal(n1::Node, n2::Node) = hash(n1)==hash(n2)


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
        dof.vals[name] = 0.0 # sets default value
    end
    return nothing
end


# Index operator for node to get a dof
# function Base.getindex(node::Node, s::Symbol)
#     return node.dofdict[s]
# end

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

"""
    $(SIGNATURES)

Get the maximum dimension of of a `nodes` collection.
"""
function getndim(nodes::Array{Node,1})
    sumz = map( node->abs(node.coord.z), nodes ) |> sum
    if sumz==0
        sumy = map( node->abs(node.coord.y), nodes ) |> sum
        ndim = sumy==0 ? 1 : 2
    else
        ndim = 3
    end

    return ndim
end


function Base.getindex(nodes::Vector{Node}, filters::NTuple; kwargs...)
    # used in nodes[ (x==0, y==1) ]
    return getindex(nodes, filters...; kwargs...)
end


# Index operator for an collection of nodes
function Base.getindex(
    nodes::Array{Node,1}, 
    filters::Union{Expr,Symbolic,String, Vector{Float64}}...;
    invert = false
    )
    
    filtered = collect(1:length(nodes))

    for filter in filters
        if typeof(filter) in (Expr, Symbolic)
            fnodes = nodes[filtered]

            T = Bool[]
            for node in fnodes
                x, y, z = node.coord.x, node.coord.y, node.coord.z
                push!(T, evaluate(filter, x=x, y=y, z=z))
            end
    
            filtered = filtered[T]
        elseif filter isa String
            filtered = [ i for i in filtered if nodes[i].tag==filter ]
        elseif filter isa Array{Float64}
            X = Vec3(filter)
            T = Bool[]
            for i in filtered
                push!(T, norm(nodes[i].coord-X) < 1e-8)
            end
            filtered = filtered[T]
        end
    end

    if invert
        filtered = setdiff(1:length(nodes), filtered)
    end

    return nodes[filtered]
end


function getnodes(nodes::Array{Node,1}, P::AbstractArray{<:Real})
    R = Node[]
    X = Vec3(P)
    for node in nodes
        if norm(X-node.coord) < 1e-8
            push!(R, node)
        end
    end
    return R
end


function getfromcoords(nodes::Array{Node,1}, P::AbstractArray{<:Real})
    R = Node[]
    X = Vec3(P)
    for node in nodes
        norm(X-node.coord) < 1e-8 && push!(R, node)
    end
    return R
end


# Get node coordinates for a collection of nodes as a matrix
function getcoords(nodes::Array{Node,1}, ndim=3)
    nnodes = length(nodes)
    return [ nodes[i].coord[j] for i in 1:nnodes, j=1:ndim]
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
    return getfield(nodes, s)
end


function nearest(nodes::Array{Node,1}, coord)
    n = length(nodes)
    D = zeros(n)
    X = Vec3(coord)

    for (i,node) in enumerate(nodes)
        D[i] = norm(X-node.coord)
    end

    return nodes[sortperm(D)[1]]
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