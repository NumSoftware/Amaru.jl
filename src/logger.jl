# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# Log structs and functions
# =========================

abstract type AbstractLogger end


# Node logger
# ===========

mutable struct NodeLogger<:AbstractLogger
    filename ::String
    filter   ::Union{Symbol,String,Expr}
    table    ::DTable
    node     ::Node

    function NodeLogger(filename::String="")
        return new(filename, :(), DTable())
    end

    #NodeLogger(node::Node, filename::String="") = new(:(), node, filename, DTable())
end

function setup_logger!(domain, filter, logger::NodeLogger)
    logger.filter = filter
    #logger.filter==:() && return
    nodes = domain.nodes[filter]
    n = length(nodes)
    n == 0 && @warn "setup_logger: No nodes found for expression:" filter=logger.filter
    n >  1 && @info "setup_logger: More than one node match expression:" filter=logger.filter
    n >= 1 && (logger.node = nodes[1])
    logger.filter = filter
    return nothing
end

function update_logger!(logger::NodeLogger, env::ModelEnv)
    isdefined(logger, :node) || return

    vals = node_vals(logger.node)
    env.transient && (vals[:t] = env.t)
    #env.transient && (vals = OrderedDict(:t=>env.t, vals...))
    push!(logger.table, vals)

    save(logger)
end


# Ip logger
# =========


mutable struct IpLogger<:AbstractLogger
    filename ::String
    filter   ::Union{Symbol,String,Expr}
    table    ::DTable
    ip       ::Ip

    function IpLogger(filename::String="")
        return new(filename, :(), DTable())
    end

    #IpLogger(ip::Ip, filename::String="") = new(:(), ip, filename, DTable())
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


function update_logger!(logger::IpLogger, env::ModelEnv)
    isdefined(logger, :ip) || return
    vals = ip_vals(logger.ip)
    env.transient && (vals[:t] = env.t)

    push!(logger.table, vals)
    save(logger)
end


# Logger for faces and edges
# ==========================


abstract type FacetLogger<:AbstractLogger
end


mutable struct FaceLogger<:FacetLogger
    filename ::String
    filter   ::Union{Symbol,String,Expr}
    table    ::DTable
    faces    ::Array{Face,1}
    nodes    ::Array{Node,1}

    function FaceLogger(filename::String="")
        return new(filename, :(), DTable())
    end
end


mutable struct EdgeLogger<:FacetLogger
    filename ::String
    filter   ::Union{Symbol,String,Expr}
    table    ::DTable
    edges    ::Array{Edge,1}
    nodes    ::Array{Node,1}

    function EdgeLogger(filename::String="")
        return new(filename, :(), DTable())
    end
end


function setup_logger!(domain, filter, logger::FaceLogger)
    logger.filter = filter
    #logger.filter==:() && return
    logger.faces = domain.faces[logger.filter]
    length(logger.faces) == 0 && @warn "setup_logger: No faces found for expression:" filter=logger.filter
    logger.nodes = logger.faces[:nodes]
end


function setup_logger!(domain, filter, logger::EdgeLogger)
    logger.filter = filter
    #logger.filter==:() && return
    logger.edges = domain.edges[logger.filter]
    length(logger.edges) == 0 && @warn "setup_logger: No edges found for expression:" filter=logger.filter
    logger.nodes = logger.edges[:nodes]
end


function update_logger!(logger::FacetLogger, env::ModelEnv)
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
    env.transient && (vals[:t] = env.t)

    push!(logger.table, vals)
    save(logger)
end


# Logger for a group of nodes
# ===========================


mutable struct NodeGroupLogger<:AbstractLogger
    filename ::String
    filter   ::Union{Symbol,String,Expr}
    book     ::DBook
    nodes    ::Array{Node,1}

    function NodeGroupLogger(filename::String="")
        return new(filename, :(), DBook(), [])
    end
end


function setup_logger!(domain, filter, logger::NodeGroupLogger)
    #logger.filter==:() && return
    logger.nodes = domain.nodes[filter]
    length(logger.nodes) == 0 && @warn "setup_logger: No nodes found for expression:" filter=logger.filter
    sort!(logger.nodes, by=n->sum(n.X))
    logger.filter = filter
end


function update_logger!(logger::NodeGroupLogger, env::ModelEnv)
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
    filename ::String
    filter   ::Union{Symbol,String,Expr}
    book     ::DBook
    ips      ::Array{Ip,1}

    function IpGroupLogger(filename::String="")
        return new(filename, :(), DBook(), [])
        #typeof(filter)<:String && ( filter=:(isequal(tag,$filter)) )
        #new(filter, Array{Ip,1}(), filename, by, DBook())
    end
    #function IpGroupLogger(ips::Array{Ip,1}, filename::String=""; by::Function=identity)
        #new(:(), ips, filename, by, DBook())
    #end
end


function setup_logger!(domain, filter, logger::IpGroupLogger)
    logger.filter = filter
    #logger.filter==:() && return
    logger.ips = domain.elems[:ips][logger.filter]
    length(logger.ips)==0 && @warn "setup_logger: No ips found for expression:" filter=logger.filter
    #logger.by != identity && sort!(logger.ips, by=logger.by)
    sort!(logger.ips, by=ip->sum(ip.X))
end


function update_logger!(logger::IpGroupLogger, env::ModelEnv)
    length(logger.ips) == 0 && return

    table = DTable()
    for ip in logger.ips
        push!(table, ip_vals(ip))
    end

    push!(logger.book, table)
    save(logger)
end


# Logger for a point
# ==================


mutable struct PointLogger<:AbstractLogger
    filename ::String
    filter   ::Array{Float64,1}
    table    ::DTable
    elem     ::Element
    R        ::Array{Float64,2}

    # logger = [  X1 => PointLogger()  ]
    # logger = [  [1,2,3] => PointLogger()  ]

    function PointLogger(filename::String="", n=20)
        return new(filename, ([0,0,0],[0,0,0]), DBook(), n, [], [])
    end
end


function setup_logger!(domain, filter, logger::PointLogger)
    logger.filter = filter
    X = filter
    # find cell and R
    mesh = domain.mesh
    cell = find_cell(X, msh.cells, msh.cellpartition, 1e-7, Cell[])
    coords = getcoords(cell)
    logger.R = inverse_map(cell.shape, coords, X)
    logger.elem = domain.elems[cell.id]
end


function update_logger!(logger::PointLogger, env::ModelEnv)
    data  = domain.point_data
    X = logger.filter
    map = [ n.id for n in logger.elem.nodes ]
    vals = OrderedDict()
    for (k,V) in data
        vals[k] = dot(V[map], N)
    end
    push!(logger.table, vals)
    save(logger)
end


# Logger for a segment
# ====================


mutable struct SegmentLogger<:AbstractLogger
    filename ::String
    filter   ::Tuple{Array{Float64,1},Array{Float64,1}}
    book     ::DBook
    n        ::Int # resolution
    elems    ::Array{Element,1}
    Rs       ::Array{Float64,2}

    # logger = [  (X1,X2) => SegmentLogger()  ]
    # logger = [  ([1,2,3],[1,2,3]) => SegmentLogger()  ]

    function SegmentLogger(filename::String="", n=20)
        return new(filename, ([0,0,0],[0,0,0]), DBook(), n, [], [])
    end
end


function setup_logger!(domain, filter, logger::SegmentLogger)
    logger.filter = filter
    X1, X2 = filter
    Xs = X1 .+ (X2.-X1)*range(0,1,length=logger.n)
    # find cell and R
    mesh = domain.mesh
    for X in Xs
        cell = find_cell(X, msh.cells, msh.cellpartition, 1e-7, Cell[])
        coords = getcoords(cell)
        R = inverse_map(cell.shape, coords, X)
        elem = domain.elems[cell.id]
        push!(logger.elems, elem)
        push!(logger.Rs, R)
    end
end

function update_logger!(logger::SegmentLogger, env::ModelEnv)
    data  = domain.point_data
    table = DTable(["s"; collect(keys(data))])
    table = DTable()
    X1, X2 = logger.filter
    Δs = norm(X2-X1)/(logger.n-1)
    s1 = 0.0
    for (elem,R) in zip(logger.elems,logger.Rs)
        s = s1 + Δs*(i-1)
        map = [ n.id for n in elem.nodes ]
        vals = [ s ]
        for (k,V) in data
            val = dot(V[map], N)
            push!(vals, val)
        end
        push!(table, vals)
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
        logger.table = DTable()
    else
        logger.book = DBook()
    end
end
