# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# Log structs and functions
# =========================

abstract type AbstractLogger end


# Node logger
# ===========

mutable struct NodeLogger<:AbstractLogger
    expr     :: Expr
    node     :: Node
    filename :: String
    table    :: DTable

    function NodeLogger(expr::Expr, filename::String="")
        this = new(expr)
        this.filename = filename
        this.table = DTable()
        return this
    end

    NodeLogger(node::Node, filename::String="") = new(:(), node, filename, DTable())
end

function setup_logger!(domain, logger::NodeLogger)
    logger.table = DTable()
    isdefined(logger, :node) && return
    nodes = domain.nodes[logger.expr]
    n = length(nodes)
    n == 0 && warn("setup_logger: No nodes found for expression: $(logger.expr)")
    n >  1 && warn("setup_logger: More than one node match expression: $(logger.expr)")
    n == 1 && (logger.node = nodes[1])
    return nothing
end

function update_logger!(logger::NodeLogger)
    isdefined(logger, :node) || return
    push!(logger.table, node_vals(logger.node))
    save(logger)
end


# Ip logger
# =========


mutable struct IpLogger<:AbstractLogger
    expr     :: Expr
    ip       :: Ip
    filename :: String
    table    :: DTable

    function IpLogger(expr::Expr, filename::String="")
        this = new(expr)
        this.filename = filename
        this.table = DTable()
        return this
    end

    IpLogger(ip::Ip, filename::String="") = new(:(), ip, filename, DTable())
end


function setup_logger!(domain, logger::IpLogger)
    logger.table = DTable()
    isdefined(logger, :ip) && return
    ips = domain.elems[:ips][logger.expr]
    n = length(ips)
    n == 0 && warn("setup_logger: No ips found for expression: $(logger.expr)")
    n >  1 && warn("setup_logger: More than one ip match expression: $(logger.expr)")
    n == 1 && (logger.ip = ips[1])
    return nothing
end


function update_logger!(logger::IpLogger)
    isdefined(logger, :ip) || return
    push!(logger.table, ip_vals(logger.ip))
    save(logger)
end


# Logger for faces and edges
# ==========================


mutable struct FacetLogger<:AbstractLogger
    expr     :: Expr
    facets   :: Array{<:Facet,1}
    nodes    :: Array{Node,1}
    filename :: String
    table    :: DTable

    function FacetLogger(applyto::Symbol, expr::Expr, filename::String="")
        this = new(expr)
        this.facets = applyto==:face ? Array{Face,1}() : Array{Edge,1}()
        this.filename = filename
        this.table = DTable()
        return this
    end

    function FacetLogger(facets::Array{<:Facet,1}, filename::String="")
        this = new()
        this.facets = facets
        this.nodes  = facets[:nodes]
        this.filename = filename
        this.table = DTable()
        return this
    end
end

function setup_logger!(domain, logger::FacetLogger)
    logger.table = DTable()
    length(logger.facets)>0 && return
    if eltype(logger.facets) == Face
        faces = domain.faces[logger.expr]
        length(faces) == 0 && warn("setup_logger: No faces found for expression: $(logger.expr)")
        logger.facets = faces
        logger.nodes = faces[:nodes]
    else
        edges = domain.edges[logger.expr]
        length(edges) == 0 && warn("setup_logger: No edges found for expression: $(logger.expr)")
        logger.facets = edges
        logger.nodes = edges[:nodes]
    end

end

function update_logger!(logger::FacetLogger)
    length(logger.facets)==0 && return

    tableU = DTable()
    tableF = DTable()
    for node in logger.nodes
        # TODO: valsF may be improved by calculating the internal forces components
        nvals = node_vals(node)
        valsU  = OrderedDict( dof.name => nvals[dof.name] for dof in node.dofs )
        valsF  = OrderedDict( dof.natname => nvals[dof.natname] for dof in node.dofs )
        push!(tableF, valsF)
        push!(tableU, valsU)
    end

    valsU = OrderedDict( key => mean(tableU[key]) for key in keys(tableU) ) # gets the average of essential values
    valsF = OrderedDict( key => sum(tableF[key])  for key in keys(tableF) ) # gets the sum for each component
    vals  = merge(valsU, valsF)

    push!(logger.table, vals)
    save(logger)
end


# Logger for a group of nodes
# ===========================


mutable struct NodeGroupLogger<:AbstractLogger
    expr     :: Expr
    nodes    :: Array{Node,1}
    filename :: String
    by       :: Function # used for sorting
    book     :: DBook

    NodeGroupLogger(expr::Expr, filename::String="", by::Function=identity) = new(expr, Array{Node,1}(), filename, by, DBook())
    NodeGroupLogger(nodes::Array{Node,1}, filename::String="", by::Function=identity) = new(:(), nodes, filename, by, DBook())
end


function setup_logger!(domain, logger::NodeGroupLogger)
    logger.book = DBook()
    length(logger.nodes) > 0 && return
    logger.nodes = domain.nodes[logger.expr]
    length(logger.nodes) == 0 && warn("setup_logger: No nodes found for expression: $(logger.expr)")
    logger.by != identity && sort!(logger.nodes, by=logger.by)
end


function update_logger!(logger::NodeGroupLogger)
    length(logger.nodes) == 0 && return

    table = DTable()
    for node in logger.nodes
        push!(table, node_vals(node))
    end
    push!(logger.book, table)
    save(logger)
end


# Logger for a group of ips
# =========================


mutable struct IpGroupLogger<:AbstractLogger
    expr     :: Expr
    ips      :: Array{Ip,1}
    filename :: String
    by       :: Function # used for sorting
    book     :: DBook

    IpGroupLogger(expr::Expr, filename::String="", by::Function=identity) = new(expr, Array{Ip,1}(), filename, by, DBook())
    IpGroupLogger(ips::Array{Ip,1}, filename::String="", by::Function=identity) = new(:(), ips, filename, by, DBook())
end


function setup_logger!(domain, logger::IpGroupLogger)
    logger.book = DBook()
    length(logger.ips)>0 && return
    logger.ips = domain.elems[:ips][logger.expr]
    length(logger.ips)==0 && warn("setup_logger: No ips found for expression: $(logger.expr)")
    logger.by != identity && sort!(logger.ips, by=logger.by)
    return nothing
end


function update_logger!(logger::IpGroupLogger)
    length(logger.ips) == 0 && return

    table = DTable()
    for ip in logger.ips
        push!(table, ip_vals(ip))
    end

    push!(logger.book, table)
    save(logger)
end


# Logger for a group of elements
# ==============================


mutable struct ElemGroupLogger<:AbstractLogger
    expr     :: Expr
    elems    :: Array{Element,1}
    nodes    :: Array{Node,1}
    ips      :: Array{Ip,1}
    filename :: String
    ipsbook  :: DBook
    nodesbook:: DBook

    function ElemGroupLogger(expr::Expr, filename::String="")
        this = new(expr)
        this.filename = filename
        this.ipsbook = DBook()
        this.nodesbook = DBook()
        return this
    end

    function ElemGroupLogger(elems::Array{<:Element,1}, filename::String="")
        this = new(:())
        this.elems = elems
        this.nodes = elems[:nodes]
        this.ips   = elems[:ips]
        this.ipsbook = DBook()
        this.nodesbook = DBook()
        return this
    end
end


function setup_logger!(domain, logger::ElemGroupLogger)
    logger.expr == :() && return
    elems = domain.elems[logger.expr]
    sort!(elems)
    n = length(elems)
    n == 0 && warn("setup_logger: No elems found for expression: $(logger.expr)")
    return nothing
end


function update_logger!(logger::ElemGroupLogger)
    isdefined(logger, :ip) || return
    push!(logger.table, ip_vals(logger.ip))
    save(logger)
end



# Function to create loggers
# ==========================


Logger(node::Node, filename::String="") = NodeLogger(node, filename)
Logger(ip::Ip, filename::String="") = IpLogger(ip, filename)
Logger(facet::Facet, filename::String="") = FacetLogger([facet], filename)
Logger(facets::Array{<:Facet,1}, filename::String="") = FacetLogger(facets, filename)
GroupLogger(nodes::Array{Node,1}, filename::String=""; by::Function=identity)=NodeGroupLogger(nodes, filename, by)
GroupLogger(ips::Array{Ip,1}, filename::String=""; by::Function=identity)=IpGroupLogger(ips, filename, by)

Logger(::Array{Node,1}, args...) = error("Logger: use GroupLogger function to log a group of nodes.")
Logger(::Array{Ip,1}, args...)   = error("Logger: use GroupLogger function to log a group of ips.")


function Logger(applyto::Symbol, target::Union{Expr, TagType}, filename="")
    available = (:node, :ip, :face, :edge)
    applyto in available || error("Logger: applyto shoud be one of $available. got :$applyto")

    typeof(target)<:TagType && (target = :(isequal(tag, $target)))

    if applyto==:node
        return NodeLogger(target, filename)
    elseif applyto==:ip
        return IpLogger(target, filename)
    else
        return FacetLogger(applyto, target, filename)
    end
end

function GroupLogger(applyto::Symbol, target::Union{Expr, TagType}, filename=""; by::Function=identity)
    available = (:node, :ip)
    applyto in available || error("GroupLogger: applyto shoud be one of $available. got :$applyto")

    typeof(target)<:TagType && (target = :(isequal(tag, $target)))

    if applyto==:node
        return NodeGroupLogger(target, filename)
    else
        return IpGroupLogger(target, filename)
    end
end


# Functions to save loggers
# =========================


function save(logger::AbstractLogger, filename::String; verbose=true)
    if isdefined(logger, :(table))
        save(logger.table, filename, verbose=verbose)
    else
        save(logger.book, filename, verbose=verbose)
    end
end


function save(logger::AbstractLogger)
    if logger.filename != ""
        save(logger, logger.filename, verbose=false)
    end
end


function reset!(logger::AbstractLogger)
    if isdefined(logger, :(table))
        logger.table = DTable()
    else
        logger.book = DBook()
    end
end
