# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

abstract type Facet end

mutable struct Face<:Facet
    shape::ShapeType
    nodes::Array{Node,1}
    ndim ::Integer
    oelem::Union{Element,Void}
    id::Int
    tag::TagType

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
    oelem::Union{Element,Void}
    id::Int
    tag::TagType

    function Edge(shape, nodes, ndim, tag="")
        this = new(shape, nodes, ndim)
        this.oelem  = nothing
        this.tag = tag
        return this
    end
end


function getindex(facets::Array{T,1}, filter_ex::Expr) where T<:Facet
    filter_ex.head in (:call, :&&, :||) || error("getindex: ill-formed condition while filterin face/edge")
    expr = fix_comparison_arrays(filter_ex)
    fun  = Functor(:(x,y,z,id,tag), expr)
   
    result = Array{T}(0)
    for facet in facets
        coords = nodes_coords(facet.nodes)
        x = coords[:,1]
        y = coords[:,2]
        z = coords[:,3]
        fun(x, y, z, facet.id, facet.tag) && push!(result, facet) 
    end

    return result
end


# Get all nodes from a collection of facets
function get_nodes(facets::Array{<:Facet,1})
    return collect( OrderedSet(node for facet in facets for node in facet.nodes) )
end

# Index operator for a collection of facets
function getindex{T<:Facet}(facets::Array{T,1}, s::Symbol)
    s == :all && return facets
    s == :nodes && return get_nodes(facets)
    error("Facet getindex: Invalid symbol $s")
end
