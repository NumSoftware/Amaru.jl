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

    function NodeLogger(expr::Union{Expr,TagType}, filename::String="")
        typeof(expr)<:TagType && ( expr=:(isequal(tag,$expr)) )
        this = new(expr)
        this.filename = filename
        this.table = DTable()
        return this
    end

    NodeLogger(node::Node, filename::String="") = new(:(), node, filename, DTable())
end

function setup_logger!(domain, logger::NodeLogger)
    logger.expr==:() && return
    nodes = domain.nodes[logger.expr]
    n = length(nodes)
    n == 0 && @warn "setup_logger: No nodes found for expression:" expr=logger.expr
    n >  1 && @info "setup_logger: More than one node match expression:" expr=logger.expr
    n >= 1 && (logger.node = nodes[1])
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

    function IpLogger(expr::Union{Expr,TagType}, filename::String="")
        typeof(expr)<:TagType && ( expr=:(isequal(tag,$expr)) )
        this = new(expr)
        this.filename = filename
        this.table = DTable()
        return this
    end

    IpLogger(ip::Ip, filename::String="") = new(:(), ip, filename, DTable())
end


function setup_logger!(domain, logger::IpLogger)
    logger.expr==:() && return
    ips = domain.elems[:ips][logger.expr]
    n = length(ips)
    n == 0 && @warn "setup_logger: No ips found for expression:" expr=logger.expr
    n >  1 && @info "setup_logger: More than one ip match expression:" expr=logger.expr
    n >= 1 && (logger.ip = ips[1])
    return nothing
end


function update_logger!(logger::IpLogger)
    isdefined(logger, :ip) || return
    push!(logger.table, ip_vals(logger.ip))
    save(logger)
end


# Logger for faces and edges
# ==========================


abstract type FacetLogger<:AbstractLogger
end


mutable struct FaceLogger<:FacetLogger
    expr     :: Expr
    faces    :: Array{Face,1}
    nodes    :: Array{Node,1}
    filename :: String
    table    :: DTable

    function FaceLogger(expr::Union{Expr,TagType}, filename::String="")
        typeof(expr)<:TagType && ( expr=:(isequal(tag,$expr)) )
        return new(expr, [], [], filename, DTable())
    end

    function FaceLogger(faces::Array{Face,1}, filename::String="")
        return new(:(), faces, faces[:nodes], filename, DTable())
    end
end


mutable struct EdgeLogger<:FacetLogger
    expr     :: Expr
    edges    :: Array{Edge,1}
    nodes    :: Array{Node,1}
    filename :: String
    table    :: DTable

    function EdgeLogger(expr::Union{Expr,TagType}, filename::String="")
        typeof(expr)<:TagType && ( expr=:(isequal(tag,$expr)) )
        return new(expr, [], [], filename, DTable())
    end

    function EdgeLogger(edges::Array{Edge,1}, filename::String="")
        return new(:(), edges, edges[:nodes], filename, DTable())
    end
end


function setup_logger!(domain, logger::FaceLogger)
    logger.expr==:() && return
    logger.faces = domain.faces[logger.expr]
    length(logger.faces) == 0 && @warn "setup_logger: No faces found for expression:" expr=logger.expr
    logger.nodes = logger.faces[:nodes]
end


function setup_logger!(domain, logger::EdgeLogger)
    logger.expr==:() && return
    logger.edges = domain.edges[logger.expr]
    length(logger.edges) == 0 && @warn "setup_logger: No edges found for expression:" expr=logger.expr
    logger.nodes = logger.edges[:nodes]
end


function update_logger!(logger::FacetLogger)
    length(logger.nodes)==0 && return

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

    function NodeGroupLogger(expr::Expr, filename::String=""; by::Function=identity) 
        typeof(expr)<:TagType && ( expr=:(isequal(tag,$expr)) )
        new(expr, Array{Node,1}(), filename, by, DBook()) 
    end
    function NodeGroupLogger(nodes::Array{Node,1}, filename::String=""; by::Function=identity)
        new(:(), nodes, filename, by, DBook()) 
    end
end


function setup_logger!(domain, logger::NodeGroupLogger)
    logger.expr==:() && return
    logger.nodes = domain.nodes[logger.expr]
    length(logger.nodes) == 0 && @warn "setup_logger: No nodes found for expression:" expr=logger.expr
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

    function IpGroupLogger(expr::Union{Expr,TagType}, filename::String=""; by::Function=identity)
        typeof(expr)<:TagType && ( expr=:(isequal(tag,$expr)) )
        new(expr, Array{Ip,1}(), filename, by, DBook())
    end
    function IpGroupLogger(ips::Array{Ip,1}, filename::String=""; by::Function=identity)
        new(:(), ips, filename, by, DBook())
    end
end


function setup_logger!(domain, logger::IpGroupLogger)
    logger.expr==:() && return
    logger.ips = domain.elems[:ips][logger.expr]
    length(logger.ips)==0 && @warn "setup_logger: No ips found for expression:" expr=logger.expr
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
