# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# Log structs and functions
# =========================

abstract type AbstractLogger end
abstract type SingleLogger<:AbstractLogger end
abstract type ComposedLogger<:AbstractLogger end


# Node logger
# ===========

mutable struct NodeLogger<:SingleLogger
    filename ::String
    filter   ::Union{Symbol,String,Expr}
    table    ::DataTable
    node     ::Node

    function NodeLogger(filename::String="")
        return new(filename, :(), DataTable())
    end
end

function setup_logger!(domain, filter, logger::NodeLogger)
    logger.filter = filter
    nodes = domain.nodes[filter]
    n = length(nodes)
    n == 0 && @warn "setup_logger: No nodes found for expression:" filter=logger.filter
    n >  1 && @info "setup_logger: More than one node match expression:" filter=logger.filter
    n >= 1 && (logger.node = nodes[1])
    logger.filter = filter
    return nothing
end

function update_logger!(logger::NodeLogger, domain)
    isdefined(logger, :node) || return

    vals = node_vals(logger.node)
    domain.env.transient && (vals[:t] = domain.env.t)
    push!(logger.table, vals)

    save(logger)
end


# Ip logger
# =========


mutable struct IpLogger<:SingleLogger
    filename ::String
    filter   ::Union{Symbol,String,Expr}
    table    ::DataTable
    ip       ::Ip

    function IpLogger(filename::String="")
        return new(filename, :(), DataTable())
    end
end


function setup_logger!(domain, filter, logger::IpLogger)
    logger.filter = filter
    #logger.filter==:() && return
    ips = domain.elems[:ips][filter]
    n = length(ips)
    n == 0 && @warn "setup_logger: No ips found for expression:" filter=logger.filter
    n >  1 && @info "setup_logger: More than one ip match expression:" filter=logger.filter
    n >= 1 && (logger.ip = ips[1])
    return nothing
end


function update_logger!(logger::IpLogger, domain)
    isdefined(logger, :ip) || return

    vals = ip_vals(logger.ip)
    domain.env.transient && (vals[:t] = domain.env.t)

    push!(logger.table, vals)
    save(logger)
end


# Logger for faces and edges
# ==========================


abstract type FacetLogger<:SingleLogger end


mutable struct FaceLogger<:FacetLogger
    filename ::String
    filter   ::Union{Symbol,String,Expr}
    table    ::DataTable
    faces    ::Array{Face,1}
    nodes    ::Array{Node,1}

    function FaceLogger(filename::String="")
        return new(filename, :(), DataTable())
    end
end


mutable struct EdgeLogger<:FacetLogger
    filename ::String
    filter   ::Union{Symbol,String,Expr}
    table    ::DataTable
    edges    ::Array{Edge,1}
    nodes    ::Array{Node,1}

    function EdgeLogger(filename::String="")
        return new(filename, :(), DataTable())
    end
end


function setup_logger!(domain, filter, logger::FaceLogger)
    logger.filter = filter
    logger.faces = domain.faces[logger.filter]
    length(logger.faces) == 0 && @warn "setup_logger: No faces found for expression:" filter=logger.filter
    logger.nodes = logger.faces[:nodes]
end


function setup_logger!(domain, filter, logger::EdgeLogger)
    logger.filter = filter
    logger.edges = domain.edges[logger.filter]
    length(logger.edges) == 0 && @warn "setup_logger: No edges found for expression:" filter=logger.filter
    logger.nodes = logger.edges[:nodes]
end


function update_logger!(logger::FacetLogger, domain)
    length(logger.nodes)==0 && return

    tableU = DataTable()
    tableF = DataTable()
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
    domain.env.transient && (vals[:t] = domain.env.t)

    push!(logger.table, vals)
    save(logger)
end


# Logger for a group of nodes
# ===========================


mutable struct NodeGroupLogger<:ComposedLogger
    filename ::String
    filter   ::Union{Symbol,String,Expr}
    book     ::DataBook
    nodes    ::Array{Node,1}

    function NodeGroupLogger(filename::String="")
        return new(filename, :(), DataBook(), [])
    end
end


function setup_logger!(domain, filter, logger::NodeGroupLogger)
    logger.nodes = domain.nodes[filter]
    length(logger.nodes) == 0 && @warn "setup_logger: No nodes found for expression:" filter=logger.filter
    sort!(logger.nodes, by=n->sum(n.X))
    logger.filter = filter
end


function update_logger!(logger::NodeGroupLogger, domain)
    length(logger.nodes) == 0 && return

    table = DataTable()
    for node in logger.nodes
        push!(table, node_vals(node))
    end
    push!(logger.book, table)
    save(logger)
end


# Logger for a group of ips
# =========================


mutable struct IpGroupLogger<:ComposedLogger
    filename ::String
    filter   ::Union{Symbol,String,Expr}
    book     ::DataBook
    ips      ::Array{Ip,1}

    function IpGroupLogger(filename::String="")
        return new(filename, :(), DataBook(), [])
    end
end


function setup_logger!(domain, filter, logger::IpGroupLogger)
    logger.filter = filter
    logger.ips = domain.elems[:ips][logger.filter]
    length(logger.ips)==0 && @warn "setup_logger: No ips found for expression:" filter=logger.filter
    sort!(logger.ips, by=ip->sum(ip.X))
end


function update_logger!(logger::IpGroupLogger, domain)
    length(logger.ips) == 0 && return

    table = DataTable()
    for ip in logger.ips
        push!(table, ip_vals(ip))
    end

    push!(logger.book, table)
    save(logger)
end


# Logger for a point
# ==================


mutable struct PointLogger<:ComposedLogger
    filename ::String
    filter   ::Array{Float64,1}
    table    ::DataTable
    elem     ::Element
    R        ::Array{Float64,1}

    # logger = [  X1 => PointLogger()  ]
    # logger = [  [1,2,3] => PointLogger()  ]

    function PointLogger(filename::String="", n=20)
        return new(filename, [0,0,0], DataTable())
    end
end


function setup_logger!(domain, filter, logger::PointLogger)
    filter isa Array{<:Number,1} || error("setup_logger!: Cannot set PointLogger. Filter should be a coordinates array.")
    logger.filter = filter
    X = filter
    # find cell and R
    mesh = domain.mesh
    cell = find_cell(X, mesh.cells, mesh._cellpartition, 1e-7, Cell[])
    cell==nothing && error("setup_logger!: Cannot set PointLogger. Coordinate ($X) outside mesh.")
    logger.R = inverse_map(cell, X)
    logger.elem = domain.elems[cell.id]
end


function update_logger!(logger::PointLogger, domain)
    data  = domain.point_data
    X = logger.filter
    N = logger.elem.shape.func(logger.R)
    map = [ n.id for n in logger.elem.nodes ]
    vals = OrderedDict()
    for (k,V) in data
        size(V,2)==1 || continue
        vals[k] = dot(V[map], N)
    end
    domain.env.transient && (vals[:t] = domain.env.t)
    push!(logger.table, vals)
    save(logger)
end


# Logger for a segment
# ====================


mutable struct SegmentLogger<:ComposedLogger
    filename ::String
    filter   ::Tuple{Array{Float64,1},Array{Float64,1}}
    book     ::DataBook
    n        ::Int # resolution
    elems    ::Array{Element,1}
    Rs       ::Array{Array{Float64,1},1}

    # logger = [  (X1,X2) => SegmentLogger()  ]
    # logger = [  ([1,2,3],[1,2,3]) => SegmentLogger()  ]

    function SegmentLogger(filename::String="", n=20)
        return new(filename, ([0,0,0],[0,0,0]), DataBook(), n, [], [])
    end
end


function setup_logger!(domain, filter, logger::SegmentLogger)
    logger.filter = filter
    X1, X2 = filter
    Δ = (X2-X1)/(logger.n-1)
    Xs = [ X1.+ Δ*i for i=0:logger.n-1]

    # find cells and Rs
    mesh = domain.mesh
    for X in Xs
        cell = find_cell(X, mesh.cells, mesh._cellpartition, 1e-7, Cell[])
        cell==nothing && error("setup_logger!: Cannot set SegmentLogger. Coordinate ($X) outside mesh.")
        coords = getcoords(cell)
        R = inverse_map(cell.shape, coords, X)
        elem = domain.elems[cell.id]
        push!(logger.elems, elem)
        push!(logger.Rs, R)
    end
end


function update_logger!(logger::SegmentLogger, domain)
    data  = domain.point_data
    labels = [ k for (k,V) in data if size(V,2)==1 ]
    table = DataTable(["s"; labels])
    X1, X2 = logger.filter
    Δs = norm(X2-X1)/(logger.n-1)
    s = 0.0
    for (elem,R) in zip(logger.elems,logger.Rs)
        N = elem.shape.func(R)
        map = [ n.id for n in elem.nodes ]
        vals = [ s ]
        for (k,V) in data
            size(V,2)==1 || continue
            val = dot(V[map], N)
            push!(vals, val)
        end
        push!(table, vals)
        s += Δs
    end

    push!(logger.book, table)
    save(logger)
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
        logger.table = DataTable()
    else
        logger.book = DataBook()
    end
end
