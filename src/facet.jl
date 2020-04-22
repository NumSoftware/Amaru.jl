# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

abstract type Facet end

mutable struct Face<:Facet
    shape::ShapeType
    nodes::Array{Node,1}
    ndim ::Integer
    oelem::Union{Element,Nothing}
    id::Int
    tag::String

    function Face(shape, nodes, ndim, tag="")
        this = new(shape, nodes, ndim)
        this.oelem  = nothing
        this.tag = tag
        return this
    end
end


mutable struct Edge<:Facet
    shape::ShapeType
    nodes::Array{Node,1}
    ndim ::Integer
    oelem::Union{Element,Nothing}
    id::Int
    tag::String

    function Edge(shape, nodes, ndim, tag="")
        this = new(shape, nodes, ndim)
        this.oelem  = nothing
        this.tag = tag
        return this
    end
end


function Base.getindex(facets::Array{T,1}, filter_ex::Expr) where T<:Facet
    nodes = facets[:nodes]
    nodemap = zeros(Int, maximum(node.id for node in nodes) )
    B = Bool[]
    for (i,node) in enumerate(nodes)
        nodemap[node.id] = i
        x, y, z = node.coord
        push!(B, eval_arith_expr(filter_ex, x=x, y=y, z=z))
    end

    R = T[]
    for facet in facets
        all( B[nodemap[node.id]] for node in facet.nodes ) && push!(R, facet)
    end
    return R
end

# Index operator for a collection of elements using a string
function Base.getindex(facets::Array{<:Facet,1}, tag::String)
    return [ facet for facet in facets if facet.tag==tag ]
end


# Get all nodes from a collection of facets
function get_nodes(facets::Array{<:Facet,1})
    return collect( OrderedSet(node for facet in facets for node in facet.nodes) )
end

# Index operator for a collection of facets
function Base.getindex(facets::Array{T,1}, s::Symbol) where T<:Facet
    s == :all && return facets
    s == :nodes && return get_nodes(facets)
    error("Facet getindex: Invalid symbol $s")
end
