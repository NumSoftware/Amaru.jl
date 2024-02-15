# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# Log structs and functions


abstract type AbstractMonitor end
@inline Base.:(<<)(a, b::AbstractMonitor) = return a=>b
# @inline Base.:(=>)(a, b::AbstractMonitor) = return (a, b)

# Ip Monitor

"""
    $(TYPEDEF)

Monitor type for a single integration point.
"""
mutable struct IpMonitor<:AbstractMonitor
    expr     ::Union{Symbol,Expr}
    val      ::Number
    filename ::String
    filter   ::Union{Int,Symbol,String,Expr,Symbolic}
    table    ::DataTable
    ip       ::Ip

    @doc """
    $(TYPEDSIGNATURES)
    
    Returns a IpMonitor for a given expression or symbol (`expr`). This monitor will
    record the corresponding data during the finite element analysis.
    The data recorded can be stored in the optional `filename`.
        

    # Example
    ```julia-repl
    julia> model = randmodel(2,2)

    julia> monitors = [
        :(x==0.5 && y==1.0) => IpMonitor(:sxx, "sxx.dat")
        :(x==0.5 && y==0.0) => IpMonitor(:eyy, "eyy.dat")
    ]

    julia> addmonitors!(model, monitors)
    ```
    """
    function IpMonitor(expr::Union{Symbol,Expr}, filename::String="")
        return new(expr, NaN, filename, :(), DataTable())
    end
end


# IpMonitor
function setup_monitor!(model, filter, monitor::IpMonitor)
    monitor.filter = filter
    if filter isa Integer
        monitor.ip = model.elems.ips[filter]
        return
    end
    ips = model.elems.ips[filter]
    n = length(ips)
    n == 0 && warn("setup_monitor!: No ips found for filter expression: $(monitor.filter)")
    n >  1 && notify("setup_monitor!: More than one ip match filter expression: $(monitor.filter)")
    n >= 1 && (monitor.ip = ips[1])
    return nothing
end


function update_monitor!(monitor::IpMonitor, model; flush=true)
    isdefined(monitor, :ip) || return success()

    state       = ip_vals(monitor.ip)
    monitor.val = eval_arith_expr(monitor.expr; state...)

    # data = OrderedDict(:val=>monitor.val)
    # data = OrderedDict(:stage=> model.env.stage, :T=>model.env.T)
    # model.env.transient && (vals[:t] = model.env.t)

    monitor.vals[:stage] = model.env.stage
    monitor.vals[:T]     = model.env.T
    model.env.transient && (monitor.vals[:t]=model.env.t)
    push!(monitor.table, data)

    if monitor.filename!="" && flush
        filename = joinpath(model.env.outdir, monitor.filename)
        save(monitor.table, filename, quiet=true)
    end
    
    return success()
end

# Node Monitor

"""
    $(TYPEDEF)

Monitor type for a single node.
"""
mutable struct NodeMonitor<:AbstractMonitor
    expr     ::Union{Symbol,Expr}
    vals     ::OrderedDict # stores current values
    filename ::String
    filter   ::Union{Int,Symbol,String,Expr,Symbolic}
    table    ::DataTable
    node     ::Node

    @doc """
    $(TYPEDSIGNATURES)
    
    Returns a NodeMonitor for a given expression or symbol (`expr`). This monitor will
    record the corresponding data during the finite element analysis.
    The data recorded can be stored in the optional `filename`.
        
    # Example
    ```julia-repl
    julia> model = randmodel(2,2)

    julia> monitors = [
        :(x==0.5 && y==1.0) => NodeMonitor(:fx, "fx.dat")
        :(x==0.5 && y==0.0) => NodeMonitor(:uy, "uy.dat")
    ]

    julia> addmonitors!(model, monitors)
    ```
    """
    function NodeMonitor(expr::Union{Symbol,Expr}, filename::String="")
        if expr isa Symbol || expr.head == :call
            expr = :($expr,)
        end

        return new(expr, OrderedDict(), filename, :(), DataTable())
    end
end


function output(monitor::AbstractMonitor)
    out = IOBuffer()
    head = "  " * string(typeof(monitor)) * " " * replace(string(monitor.filter), " " => "")
    print(out, head, "\e[K")

    first = true
    for (k,v) in monitor.vals
        k in (:stage, :T) && continue
        first || print(out, " "^length(head))
        first = false

        # v = v isa Real ? round(v, sigdigits=5) : v
        if v isa Real && !(v isa Bool)
            v = round(v, sigdigits=5)
        end
        println(out, "   ", replace(string(k), " " => ""), " : ", v, "\e[K")
    end
    str = String(take!(out))
    # return replace(str, "\n" =>)
    return str
end


function setup_monitor!(model, filter, monitor::NodeMonitor)
    monitor.filter = filter
    if filter isa Integer
        monitor.node = model.elems.nodes[filter]
        return
    end
    nodes = model.elems.nodes[filter]
    n = length(nodes)
    n == 0 && warn("setup_monitor!: No nodes found for filter expression: $(monitor.filter)")
    n >  1 && notify("setup_monitor!: More than one node match filter expression: $(monitor.filter)")
    n >= 1 && (monitor.node = nodes[1])
    return nothing
end


function update_monitor!(monitor::NodeMonitor, model; flush=true)
    isdefined(monitor, :node) || return success()

    for expr in monitor.expr.args
        state = node_vals(monitor.node) 
        monitor.vals[expr] = eval_arith_expr(expr; state...)
    end

    monitor.vals[:stage] = model.env.stage
    monitor.vals[:T]     = model.env.T
    model.env.transient && (monitor.vals[:t]=model.env.t)

    push!(monitor.table, monitor.vals)

    if monitor.filename!="" && flush
        filename = joinpath(model.env.outdir, monitor.filename)
        save(monitor.table, filename, quiet=true)
    end

    return success()
end


# Monitor for a group of ips


mutable struct IpGroupMonitor<:AbstractMonitor
    expr    ::Expr
    vals    ::OrderedDict
    filename::String
    filter  ::Union{Symbol,String,Expr,Symbolic}
    table   ::DataTable
    ips     ::Array{Ip,1}

    function IpGroupMonitor(expr::Expr, filename::String="")
        if expr.head == :call
            expr = :($expr,)
        end

        return new(expr, OrderedDict(), filename, :(), DataTable(), [])
    end
end


function setup_monitor!(model, filter, monitor::IpGroupMonitor)
    monitor.filter = filter
    monitor.ips = model.elems.ips[monitor.filter]
    length(monitor.ips)==0 && warn("setup_monitor!: No ips found for filter expression: ", monitor.filter)
end


function update_monitor!(monitor::IpGroupMonitor, model; flush=true)
    length(monitor.ips) == 0 && return success()

    for expr in monitor.expr.args
        vals = []
        for ip in monitor.ips
            state = ip_vals(ip) 
            # @s expr
            # @s state
            val = eval_arith_expr(expr; state...)
            push!(vals, val)
        end
        if expr isa Symbol
            lo, hi = round.(extrema(vals), sigdigits=5)
            monitor.vals[:(min($expr))] = lo
            monitor.vals[:(max($expr))] = hi
        else
            val = any(vals)
            monitor.vals[expr] = val
        end
    end

    monitor.vals[:stage] = model.env.stage
    monitor.vals[:T]     = model.env.T
    model.env.transient && (monitor.vals[:t]=model.env.t)

    push!(monitor.table, monitor.vals)

    if monitor.filename!="" && flush
        filename = joinpath(model.env.outdir, monitor.filename)
        save(monitor.table, filename, quiet=true)
    end

    return success()
end


# Logger for the sum of nodes


"""
    $(TYPEDEF)

Monitor type tha represent the resultant of a group of nodes.
"""
mutable struct NodeSumMonitor<:AbstractMonitor
    expr     ::Union{Symbol,Expr} # watches
    vals     ::OrderedDict # current values
    filename ::String
    filter   ::Union{Symbol,String,Expr,Symbolic}
    table    ::DataTable
    nodes    ::Array{Node,1}
    stopexpr ::Expr
    # _extra   ::Expr # extra symbols :(fx_min, fx_max, stress=fx/10), etc

    @doc """
    $(TYPEDSIGNATURES)
    
    Returns a NodeSumMonitor for a given expression or symbol (`expr`). This monitor will
    record the corresponding data during the finite element analysis.
    The data recorded can be stored in the optional `filename`.
    The `stop` expresion can be used to stop the analysis if the expresion becomes true.

    # Example
    ```julia-repl
    julia> model = randmodel(2,2)

    julia> monitors = [
        :(y==1.0) => NodeSumMonitor(:fx, "fx.dat")
        :(y==0.0) => NodeSumMonitor(:fy, "fy.dat"; stop = :(uy>0.01))
    ]

    julia> addmonitors!(model, monitors)
    ```
    """
    function NodeSumMonitor(expr::Union{Symbol,Expr}, filename::String=""; stop=Expr=:())
        if expr isa Symbol || expr.head == :call || expr.head == :(=)
            expr = :($expr,)
        end
        if stop.head == :call
            stop = :($stop,)
        end
        # add _min and _max extra values
        # extra = :()
        # for ex in [ expr.args; stop.args ]
        #     for var in getvars(ex)
        #         varstr = string(var)
        #         length(varstr)>4 || continue
        #         _, optstr = split(varstr, '_')
        #         optstr in ("min","max") || continue
        #         extra.args = union(extra.args, [var]) # add extra var such as fx_min
        #         expr.args  = union(expr.args, [var])  # add extra var such as fx
        #     end
        # end

        return new(expr, OrderedDict(), filename, :(), DataTable(), Node[], stop)
    end
end


function setup_monitor!(model, filter, monitor::NodeSumMonitor)
    monitor.filter = filter
    monitor.nodes = model.nodes[filter]
    length(monitor.nodes) == 0 && warn("setup_monitor: No nodes found for filter expression: ", monitor.filter)
end


function update_monitor!(monitor::NodeSumMonitor, model; flush=true)
    length(monitor.nodes) == 0 && return success()

    tableU = DataTable()
    tableF = DataTable()
    for node in monitor.nodes
        nvals = node_vals(node)
        valsU  = OrderedDict( dof.name => nvals[dof.name] for dof in node.dofs )
        valsF  = OrderedDict( dof.natname => nvals[dof.natname] for dof in node.dofs )
        push!(tableF, valsF)
        push!(tableU, valsU)
    end

    valsU = OrderedDict( Symbol(key) => mean(tableU[key]) for key in keys(tableU) ) # gets the average of essential values
    valsF = OrderedDict( Symbol(key) => sum(tableF[key])  for key in keys(tableF) ) # gets the sum for each component
    state = merge(valsU, valsF)

    # update state with definitinos
    for ex in monitor.expr.args
        if ex isa Expr && ex.head == :(=)
            var = ex.args[1]
            state[var] = eval_arith_expr(ex.args[2]; state...)
        end
    end

    # update state with _max and _min
    vars = union(getvars(monitor.expr), getvars(monitor.stopexpr))
    extravals = OrderedDict()
    for var in vars
        contains(string(var), r"_max|_min") || continue
        keystr, optstr = split(string(var), '_')
        key = Symbol(keystr)

        if size(monitor.table)[1] == 0
            state[var] = state[key]
        else
            fun = optstr=="min" ? min : max
            state[var] = fun(state[key], monitor.table[var][end])
        end
        extravals[var] = state[var]
    end

    # eval expressions
    for expr in monitor.expr.args
        if expr isa Expr && expr.head == :(=)
            expr = expr.args[1]
        end
        monitor.vals[expr] = eval_arith_expr(expr; state...)
    end

    merge!(monitor.vals, extravals)
    monitor.vals[:stage] = model.env.stage
    monitor.vals[:T]     = model.env.T
    model.env.transient && (monitor.vals[:t]=model.env.t)
    push!(monitor.table, monitor.vals)

    if monitor.filename!="" && flush
        filename = joinpath(model.env.outdir, monitor.filename)
        save(monitor.table, filename, quiet=true)
    end

    # eval stop expressions
    for expr in monitor.stopexpr.args

        # @show expr
        # @show eval_arith_expr(expr; state...)
        # @show eval_arith_expr(:fx; state...)
        # @show eval_arith_expr(:fx_max; state...)
        if eval_arith_expr(expr; state...)
            # @show "STOPPPPPPPP"
            return failure("Stop condition at NodeSumMonitor ($expr)")
        end
        # error()
    end

    return success()
end
