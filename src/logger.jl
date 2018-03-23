# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

abstract type AbstractLogger end

mutable struct Logger<:AbstractLogger
    applyto  :: Symbol
    expr     :: Expr
    objects  :: Array{T,1} where T <: Union{Node, Ip}
    filename :: String
    table    :: DTable

    function Logger(applyto::Symbol, target::Union{Expr, TagType}, filename="")
        @assert(applyto in (:node, :ip, :faces, :edges) )
        this = new()
        this.applyto  = applyto
        this.table    = DTable()
        #this.objects  = []
        this.filename = filename

        ty = typeof(target)
        if ty==Expr
            this.expr = target
        #elseif ty==Int || ty==Array{Int,1}
            #this.expr = :(id in $target)
        elseif ty<:Int || ty<:String
            this.expr = :(isequal(tag, $target))
        end
        return this
    end
end

mutable struct GroupLogger<:AbstractLogger
    applyto  :: Symbol
    expr     :: Expr
    objects  :: Array{T,1} where T <: Union{Node, Ip}
    filename :: String
    by       :: Function # used for sorting
    book     :: DBook

    function GroupLogger(applyto::Symbol, target::Union{Expr, TagType}, filename=""; by::Function=identity)
        @assert(applyto in (:nodes, :ips, :node, :ip) )
        this = new()
        this.applyto  = applyto
        this.book     = DBook()
        #this.objects  = []
        this.filename = filename
        this.by       = by

        ty = typeof(target)
        if ty==Expr
            this.expr = target
        elseif ty==Int || ty==Array{Int,1}
            this.expr = :(id in $target)
        elseif ty==String
            this.expr = :(isequal(tag, $target))
        end
        return this
    end
end


###############################################################################
# Functions to setup loggers

function set_logger(domain, logger::Logger)
    if logger.applyto == :node
        nodes = domain.nodes[logger.expr]
        logger.objects = nodes
        if length(nodes) == 0
            warn("set_logger: No nodes found for expression: $(logger.expr)")
        elseif length(nodes) > 1
            warn("set_logger: More than one node match expression: $(logger.expr)")
        end
    elseif logger.applyto == :ip
        ips = domain.elems[:ips][logger.expr]
        logger.objects = ips
        if length(ips) == 0
            warn("set_logger: No ips found for expression: $(logger.expr)")
        elseif length(ips) > 1
            warn("set_logger: More than one ip match expression: $(logger.expr)")
        end
    elseif logger.applyto == :faces
        faces = domain.faces[logger.expr]
        length(faces) == 0 && warn("set_logger: No faces found for expression: $(logger.expr)")
        logger.objects = faces[:nodes]
    elseif logger.applyto == :edges
        edges = domain.edges[logger.expr]
        length(edges) == 0 && warn("set_logger: No edges found for expression: $(logger.expr)")
        logger.objects = edges[:nodes]
    end

end


function set_logger(domain, logger::GroupLogger)
    if logger.applyto == :nodes
        logger.objects = domain.nodes[logger.expr]
        length(logger.objects) == 0 && warn("set_logger: No nodes found for expression: $(logger.expr)")

        if logger.by != identity
            logger.objects = sort(logger.objects, by=logger.by)
        end
    elseif logger.applyto == :ips
        logger.objects = domain.elems[:ips][logger.expr]
        length(logger.objects) == 0 && warn("set_logger: No ips found for expression: $(logger.expr)")

        if logger.by != identity
            logger.objects = sort(logger.objects, by=logger.by)
        end
    end
end


###############################################################################
# Functions to update loggers

function save(logger::AbstractLogger, filename::String; verbose=true)
    if isdefined(logger, :(table))
        save(logger.table, filename, verbose=verbose)
    else
        save(logger.book, filename, verbose=verbose)
    end
end

function save(logger::AbstractLogger)
    if logger.filename != ""
        try   save(logger, logger.filename, verbose=false)
        catch warn("Problem writing file ", logger.filename)
        end
    end
end


function update_logger!(logger::Logger)
    length(logger.objects)==0 && return
    local vals
    if logger.applyto == :node
        node = logger.objects[1]
        vals = node_vals(node)
    elseif logger.applyto == :ip
        ip   = logger.objects[1]
        vals = ip_state_vals(ip.owner.mat, ip.data)
    elseif logger.applyto in (:faces, :edges)
        tableU = DTable()
        tableF = DTable()
        for node in logger.objects
            # TODO: valsF may be improved by calculating the internal forces components
            nvals = node_vals(node)
            valsU  = Dict( dof.name => nvals[dof.name] for dof in node.dofs )
            valsF  = Dict( dof.natname => nvals[dof.natname] for dof in node.dofs )
            push!(tableF, valsF)
            push!(tableU, valsU)
        end

        valsU = Dict( key => mean(tableU[key]) for key in keys(tableU) ) # gets the average of essential values
        valsF = Dict( key => sum(tableF[key])  for key in keys(tableF) ) # gets the sum for each component
        vals  = merge(valsU, valsF)
    end
    push!(logger.table, vals)
    save(logger)
end


function update_logger!(logger::GroupLogger)
    length(logger.objects)==0 && return
    table = DTable()

    if logger.applyto == :nodes
        for node in logger.objects
            vals = node_vals(node)
            push!(table, vals)
        end
    elseif logger.applyto == :ips
        for ip in logger.objects
            vals = ip_vals(ip)  # includes ip global coordinates
            push!(table, vals)
        end
    end
    push!(logger.book, table)
    save(logger)
end

function update_loggers!(domain)
    for logger in domain.loggers
        update_logger!(logger)
    end
end



