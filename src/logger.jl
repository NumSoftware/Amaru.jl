# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# Log structs and functions
# =========================

abstract type AbstractLogger end
abstract type SingleLogger<:AbstractLogger end
abstract type MultiLogger<:AbstractLogger end

@inline Base.:(<<)(a, b::AbstractLogger) = return (a, b)
@inline Base.:(=>)(a, b::AbstractLogger) = return (a, b)

# Node logger
# ===========

mutable struct NodeLogger<:SingleLogger
    filename ::String
    filter   ::Any
    table    ::DataTable
    node     ::Node

    function NodeLogger(filename::String="")
        return new(filename, :(), DataTable())
    end
end


function setup_logger!(ana::Analysis, filter, logger::NodeLogger)
    logger.filter = filter

    if filter isa AbstractArray
        X = Vec3(filter)
        x, y, z = X
        nodes = ana.model.nodes[:(x==$x && y==$y && z==$z)]
        n = length(nodes)

        if n==0
            logger.node = nearest(ana.model.nodes, X)
            notify("setup_logger: No node found at $(logger.filter). Picking the nearest at $(logger.node.coord)")
        else
            logger.node = nodes[1]
        end
        logger.filename = getfullpath(ana.sctx.outdir, logger.filename)

        return
    end

    nodes = ana.model.nodes[filter]
    n = length(nodes)
    n == 0 && warn("setup_logger: No nodes found for filter expression: ", logger.filter)
    n >  1 && notify("setup_logger: More than one node match filter expression: ", logger.filter)
    n >= 1 && (logger.node = nodes[1])
    logger.filename = getfullpath(ana.sctx.outdir, logger.filename)

    return nothing
end


function update_logger!(logger::NodeLogger, ana::Analysis)
    isdefined(logger, :node) || return

    vals = node_vals(logger.node)
    ana.sctx.transient && (vals[:t] = ana.sctx.t)
    push!(logger.table, vals)
end


# Ip logger
# =========

mutable struct IpLogger<:SingleLogger
    filename ::String
    filter   ::Any
    table    ::DataTable
    ip       ::Ip

    function IpLogger(filename::String="")
        return new(filename, :(), DataTable())
    end
end


function setup_logger!(ana::Analysis, filter, logger::IpLogger)
    logger.filter = filter
    if filter isa Integer
        logger.ip = ana.model.elems.ips[filter]
        return
    end
    if filter isa AbstractArray
        X = Vec3(filter)
        x, y, z = X
        ips = ana.model.elems.ips[:(x==$x && y==$y && z==$z)]
        n = length(ips)

        if n==0
            logger.ip = nearest(ana.model.elems.ips, X)
            X = round.(logger.ip.coord, sigdigits=5)
            notify("setup_logger: No ip found at $(logger.filter). Picking the nearest at $X")
        else
            logger.ip = ips[1]
        end
        return
    end

    ips = ana.model.elems.ips[filter]
    n = length(ips)
    n == 0 && warn("setup_logger: No ips found for filter expression: $(logger.filter)")
    n >  1 && notify("setup_logger: More than one ip match filter expression: $(logger.filter)")
    n >= 1 && (logger.ip = ips[1])

    logger.filename = getfullpath(ana.sctx.outdir, logger.filename)

    return nothing
end


function update_logger!(logger::IpLogger, ana::Analysis)
    isdefined(logger, :ip) || return

    vals = ip_vals(logger.ip)
    vals[:out] = ana.sctx.out
    vals[:T] = ana.sctx.T
    ana.sctx.transient && (vals[:t] = ana.sctx.t)

    push!(logger.table, vals)

    # if logger.filename!="" && flush
    #     filename = joinpath(ana.sctx.outdir, logger.filename)
    #     save(logger, filename, quiet=true)
    # end
end


# Logger for faces and edges



abstract type FacetLogger<:SingleLogger end


mutable struct FaceLogger<:FacetLogger
    filename ::String
    filter   ::Any
    table    ::DataTable
    faces    ::Array{CellFace,1}
    nodes    ::Array{Node,1}

    function FaceLogger(filename::String="")
        return new(filename, :(), DataTable())
    end
end


mutable struct EdgeLogger<:FacetLogger
    filename ::String
    filter   ::Any
    table    ::DataTable
    edges    ::Array{CellEdge,1}
    nodes    ::Array{Node,1}

    function EdgeLogger(filename::String="")
        return new(filename, :(), DataTable())
    end
end

FacesSumLogger = FaceLogger
SurfaceLogger = FaceLogger
EdgesSumLogger = EdgeLogger

export FacesSumLogger, EdgesSumLogger, NodesSumLogger


function setup_logger!(ana::Analysis, filter, logger::FaceLogger)
    if length(ana.model.elems.joints) > 0
        warn("setup_logger: Using FaceLogger in a mesh with joint elements may provide wrong results.")
    end

    logger.filter = filter
    logger.faces = ana.model.faces[logger.filter]
    length(logger.faces) == 0 && warn("setup_logger: No faces found for filter expression: ", logger.filter)
    logger.nodes = logger.faces.nodes
    logger.filename = getfullpath(ana.sctx.outdir, logger.filename)

    return nothing
end


function setup_logger!(ana::Analysis, filter, logger::EdgeLogger)
    logger.filter = filter
    logger.edges = ana.model.edges[logger.filter]
    length(logger.edges) == 0 && warn("setup_logger: No edges found for filter expression: ", logger.filter)
    logger.nodes = logger.edges.nodes
    logger.filename = getfullpath(ana.sctx.outdir, logger.filename)

    return nothing
end


function update_logger!(logger::FacetLogger, ana::Analysis)
    length(logger.nodes)==0 && return

    tableU = DataTable()
    tableF = DataTable()
    for node in logger.nodes
        # TODO: valsF may be improved by calculating the internal forces components
        nvals = node_vals(node)
        valsU = OrderedDict( dof.name => nvals[dof.name] for dof in node.dofs )
        valsF = OrderedDict( dof.natname => nvals[dof.natname] for dof in node.dofs )
        push!(tableF, valsF)
        push!(tableU, valsU)
    end

    valsU = OrderedDict( Symbol(key) => mean(tableU[key]) for key in keys(tableU) ) # gets the average of essential values
    valsF = OrderedDict( Symbol(key) => sum(tableF[key])  for key in keys(tableF) ) # gets the sum for each component
    vals  = merge(valsU, valsF)
    vals[:out] = ana.sctx.out
    ana.sctx.transient && (vals[:t] = ana.sctx.t)

    push!(logger.table, vals)
end


# Logger for a pack of nodes



mutable struct NodeSumLogger<:SingleLogger
    filename ::String
    filter   ::Any
    table    ::DataTable
    nodes    ::Array{Node,1}

    function NodeSumLogger(filename::String="")
        return new(filename, :(), DataTable(), [])
    end
end

NodesSumLogger = NodeSumLogger


function setup_logger!(ana::Analysis, filter, logger::NodeSumLogger)
    logger.filter = filter
    logger.nodes = ana.model.nodes[filter]
    length(logger.nodes) == 0 && warn("setup_logger: No nodes found for filter expression: ", logger.filter)
    logger.filename = getfullpath(ana.sctx.outdir, logger.filename)

    return nothing
end


function update_logger!(logger::NodeSumLogger, ana::Analysis)
    length(logger.nodes) == 0 && return

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

    valsU = OrderedDict( Symbol(key) => mean(tableU[key]) for key in keys(tableU) ) # gets the average of essential values
    valsF = OrderedDict( Symbol(key) => sum(tableF[key])  for key in keys(tableF) ) # gets the sum for each component
    vals  = merge(valsU, valsF)
    vals[:out] = ana.sctx.out
    ana.sctx.transient && (vals[:t] = ana.sctx.t)

    push!(logger.table, vals)
end


# Logger for a group of nodes


mutable struct NodeGroupLogger<:MultiLogger
    filename ::String
    filter   ::Any
    book     ::DataBook
    nodes    ::Array{Node,1}

    function NodeGroupLogger(filename::String="")
        return new(filename, :(), DataBook(), [])
    end
end


function setup_logger!(ana::Analysis, filter, logger::NodeGroupLogger)
    logger.filter = filter
    logger.nodes  = ana.model.nodes[filter]
    length(logger.nodes) == 0 && warn("setup_logger: No nodes found for filter expression: ", logger.filter)

    # first try to sort nodes according to its elements (useful for nodes with the same coordintates)
    dists = zeros(length(logger.nodes))
    for elem in ana.model.elems
        elem.shape.family in (JOINTCELL, LINEJOINTCELL, TIPJOINTCELL) && continue
        for (i,node) in enumerate(logger.nodes)
            if node in elem.nodes
                dists[i] = sum(elem.coords)
            end
        end
    end

    # check if distances vector is valid
    if !(0.0 in dists)
        idxs = sortperm(dists)
        logger.nodes = logger.nodes[idxs]
    end

    # sort nodes
    sort!(logger.nodes, by=n->sum(n.coord))
    
    logger.filename = getfullpath(ana.sctx.outdir, logger.filename)
end


function update_logger!(logger::NodeGroupLogger, ana::Analysis)
    length(logger.nodes) == 0 && return

    table = DataTable()

    # add data
    ids = [ node.id for node in logger.nodes ]
    for (k,V) in ana.model.node_data
        table[k] = V[ids]
    end

    # add coordinates
    for (i,k) in enumerate((:x, :y, :z))
        table[k] = [ node.coord[i] for node in logger.nodes ]
    end

    # for node in logger.nodes
    #     # push!(table, node_vals(node))
    # end
    push!(logger.book, table)

    # if logger.filename!="" && flush
    #     filename = joinpath(ana.sctx.outdir, logger.filename)
    #     save(logger, filename, quiet=true)
    # end
end


# Logger for a group of ips

mutable struct IpGroupLogger<:MultiLogger
    filename ::String
    filter   ::Any
    book     ::DataBook
    ips      ::Array{Ip,1}

    function IpGroupLogger(filename::String="")
        return new(filename, :(), DataBook(), [])
    end
end


function setup_logger!(ana::Analysis, filter, logger::IpGroupLogger)
    logger.filter = filter
    logger.ips = ana.model.elems.ips[logger.filter]
    length(logger.ips)==0 && warn("setup_logger: No ips found for filter expression: ", logger.filter)
    sort!(logger.ips, by=ip->sum(ip.coord))

    logger.filename = getfullpath(ana.sctx.outdir, logger.filename)
    return nothing
end


function update_logger!(logger::IpGroupLogger, ana::Analysis)
    length(logger.ips) == 0 && return

    table = DataTable()
    for ip in logger.ips
        push!(table, ip_vals(ip))
    end

    push!(logger.book, table)
end


# Logger for a point
# ==================


mutable struct PointLogger<:SingleLogger
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


function setup_logger!(ana::Analysis, filter, logger::PointLogger)
    filter isa Array{<:Number,1} || error("setup_logger!: Cannot set PointLogger. Filter should be a coordinates array.")
    logger.filter = filter
    X = filter
    # find elem and R
    elem = find_elem(X, ana.model.elems, ana.model._elempartition)
    elem===nothing && error("setup_logger!: Cannot set PointLogger. Coordinate ($X) outside mesh.")
    logger.elem = elem
    logger.R = inverse_map(elem, X)

    logger.filename = getfullpath(ana.sctx.outdir, logger.filename)

    return nothing
end


function update_logger!(logger::PointLogger, ana::Analysis)
    data  = ana.model.node_data
    N = logger.elem.shape.func(logger.R)
    map = [ n.id for n in logger.elem.nodes ]
    vals = OrderedDict()
    for (k,V) in data
        size(V,2)==1 || continue
        vals[k] = dot(V[map], N)
    end
    vals[:out] = ana.sctx.out
    ana.sctx.transient && (vals[:t] = ana.sctx.t)
    push!(logger.table, vals)
end


# Logger for a segment
# ====================


mutable struct SegmentLogger<:MultiLogger
    filename ::String
    filter   ::Array{Float64,2}
    book     ::DataBook
    n        ::Int # resolution
    elems    ::Array{Element,1}
    Rs       ::Array{Array{Float64,1},1}

    # logger = [  (X1,X2) => SegmentLogger()  ]
    # logger = [  ([1,2,3],[1,2,3]) => SegmentLogger()  ]

    function SegmentLogger(filename::String=""; n=20)
        return new(filename, [0.0 0 0; 0 0 0], DataBook(), n, [], [])
    end
end


function setup_logger!(ana::Analysis, filter, logger::SegmentLogger)
    logger.filter = filter
    X1 = filter[1,:]
    X2 = filter[2,:]
    Δ = (X2-X1)/(logger.n-1)
    Xs = [ X1.+ Δ*i for i in 0:logger.n-1]

    # find cells and Rs
    for X in Xs
        elem = find_elem(X, ana.model.elems, ana.model._elempartition)
        elem===nothing && error("setup_logger!: Cannot set SegmentLogger. Coordinate ($X) outside mesh.")
        R = inverse_map(elem, X)
        push!(logger.elems, elem)
        push!(logger.Rs, R)
    end

    logger.filename = getfullpath(ana.sctx.outdir, logger.filename)

    return nothing
end


function update_logger!(logger::SegmentLogger, ana::Analysis)
    ndim = ana.model.ctx.ndim
    data  = ana.model.node_data
    coord_labels = ["x", "y", "z"][1:ndim]
    labels = [ k for (k,V) in data if size(V,2)==1 ]
    table = DataTable(["s"; coord_labels; labels])
    X1, X2 = logger.filter
    Δs = norm(X2-X1)/(logger.n-1)
    s = 0.0
    for (elem,R) in zip(logger.elems,logger.Rs)
        N = elem.shape.func(R)
        X = getcoords(elem)'*N
        map = [ n.id for n in elem.nodes ]
        vals = [ s; X ]
        for (k,V) in data
            size(V,2)==1 || continue
            val = dot(V[map], N)
            push!(vals, val)
        end
        push!(table, vals)
        s += Δs
    end

    push!(logger.book, table)
end


# Functions to save loggers

function save(logger::SingleLogger, filename::String; quiet=false)
    save(logger.table, filename, quiet=quiet)
end
function save(logger::MultiLogger, filename::String; quiet=false)
    save(logger.book, filename, quiet=quiet)
end
