# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# Monitor structs and functions
# =============================

abstract type AbstractMonitor end
@inline Base.:(<<)(a, b::AbstractMonitor) = return a=>b

# Ip Monitor
# ==========

"""
    $(TYPEDEF)

Monitor type for a single integration point.
"""
mutable struct IpMonitor<:AbstractMonitor
    expr     ::Union{Symbol,Expr}
    vals     ::OrderedDict # stores current values
    filename ::String
    filter   ::Union{Int,Symbol,String,Expr,Symbolic,AbstractArray}
    label    ::String
    table    ::DataTable
    ip       ::Union{Ip,Nothing}

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
        if expr isa Expr
            expr = round_floats!(expr)
        end
        if expr isa Symbol || expr.head == :call || expr.head == :(=)
            expr = :($expr,)
        end
        return new(expr, OrderedDict(), filename, :(), "", DataTable())
    end
end


function setup_monitor!(model, allitems, filter, monitor::AbstractMonitor, field::Symbol)
    item_type  = eltype(allitems)
    item_name  = lowercase(string(item_type))
    singleitem = string(field)[end] != 's' # if field ends with 's' it is a collection

    if filter isa AbstractArray
        items = getfromcoords(allitems, filter)
        n = length(items)

        if n==0
            N = nearest(allitems, filter)
            if N===nothing
                items = item_type[]
            else
                notify("setup_monitor: No ip found at $(monitor.filter). Picking the nearest at $(N.coord)")
                items = [ N ]
            end
        end
    elseif filter isa Integer
        if 1<=filter<=length(allitems)
            items = [ allitems[filter] ]
        else
            items = item_type[]
        end
    else
        items = allitems[filter]
    end

    n = length(items)
    if singleitem
        if n==0
            warn("setup_monitor!: No $(item_name)s found for filter expression: $(monitor.filter)")
            setfield!(monitor, field, nothing)
        else
            n>1 && notify("setup_monitor!: More than one $item_name match filter expression: $(monitor.filter)")
            setfield!(monitor, field, items[1])
        end
    else
        setfield!(monitor, field, items)
    end

    if filter isa Expr || filter isa Symbolic
        monitor.label = replace( string(round_floats!(getexpr(filter))), " " => "")
    elseif filter isa Integer
        monitor.label = "$item_name $filter"
    else
        monitor.label = string(filter)
    end

    monitor.filename = getfullpath(model.env.outdir, monitor.filename)

    return nothing
end


function setup_monitor!(model, filter, monitor::IpMonitor)
    setup_monitor!(model, model.elems.ips, filter, monitor, :ip)
end


function update_monitor!(monitor::IpMonitor, model)
    monitor.ip === nothing && return success()

    for expr in monitor.expr.args
        state = ip_vals(monitor.ip) 
        monitor.vals[expr] = evaluate(expr; state...)
    end

    monitor.vals[:stage] = model.env.stage
    monitor.vals[:T]     = model.env.T
    model.env.transient && (monitor.vals[:t]=model.env.t)
    push!(monitor.table, monitor.vals)

    return success()
end


function output(monitor::AbstractMonitor)
    out = IOBuffer()
    head = "  " * string(typeof(monitor)) * " " * monitor.label
    print(out, head, "\e[K")

    first = true
    for (k,v) in monitor.vals
        k in (:stage, :T, :t) && continue
        first || print(out, " "^length(head))
        first = false

        if v isa Real && !(v isa Bool)
            v = round(v, sigdigits=5)
        end
        println(out, "   ", replace(string(k), " " => ""), " : ", v, "\e[K")
    end
    str = String(take!(out))
    return str
end


# Node Monitor
# ============

"""
    $(TYPEDEF)

Monitor type for a single node.
"""
mutable struct NodeMonitor<:AbstractMonitor
    expr     ::Union{Symbol,Expr}
    vals     ::OrderedDict # stores current values
    filename ::String
    filter   ::Union{Int,Symbol,String,Expr,Symbolic,AbstractArray}
    label    ::String
    table    ::DataTable
    node     ::Union{Node,Nothing}

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

        return new(expr, OrderedDict(), filename, :(), "", DataTable())
    end
end


function setup_monitor!(model, filter, monitor::NodeMonitor)
    setup_monitor!(model, model.nodes, filter, monitor, :node)
end


function update_monitor!(monitor::NodeMonitor, model)
    monitor.node === nothing && return success()

    for expr in monitor.expr.args
        state = node_vals(monitor.node) 
        monitor.vals[expr] = evaluate(expr; state...)
    end

    monitor.vals[:stage] = model.env.stage
    monitor.vals[:T]     = model.env.T
    model.env.transient && (monitor.vals[:t]=model.env.t)

    push!(monitor.table, monitor.vals)

    return success()
end


# Monitor for a group of ips
# ==========================

mutable struct IpGroupMonitor<:AbstractMonitor
    expr    ::Expr
    vals    ::OrderedDict
    filename::String
    filter  ::Union{Symbol,String,Expr,Symbolic}
    label   ::String
    table   ::DataTable
    ips     ::Array{Ip,1}

    function IpGroupMonitor(expr::Union{Symbol,Expr}, filename::String="")
        if expr isa Expr
            expr = round_floats!(expr)
        end
        if expr isa Symbol || expr.head == :call || expr.head == :(=)
            expr = :($expr,)
        end

        return new(expr, OrderedDict(), filename, :(), "", DataTable(), [])
    end
end


function setup_monitor!(model, filter, monitor::IpGroupMonitor)
    setup_monitor!(model, model.elems.ips, filter, monitor, :ip)
end


function update_monitor!(monitor::IpGroupMonitor, model)
    length(monitor.ips) == 0 && return success()

    for expr in monitor.expr.args
        vals = []
        for ip in monitor.ips
            state = ip_vals(ip) 
            val = evaluate(expr; state...)
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

    return success()
end


# Monitor for the sum of nodes
# ============================

"""
    $(TYPEDEF)

Monitor type tha represent the resultant of a group of nodes.
"""
mutable struct NodeSumMonitor<:AbstractMonitor
    expr     ::Union{Symbol,Expr} # watches
    vals     ::OrderedDict # current values
    filename ::String
    filter   ::Union{Symbol,String,Expr,Symbolic}
    label    ::String
    table    ::DataTable
    nodes    ::Array{Node,1}
    stopexpr ::Expr

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
        if expr isa Expr
            expr = fix_expr_maximum_minimum!(round_floats!(expr))
        end
        stop = fix_expr_maximum_minimum!(stop)
        if expr isa Symbol || expr.head == :call || expr.head == :(=)
            expr = :($expr,)
        end
        if stop.head == :call
            stop = :($stop,)
        end

        return new(expr, OrderedDict(), filename, :(), "", DataTable(), Node[], stop)
    end
end


function setup_monitor!(model, filter, monitor::NodeSumMonitor)
    setup_monitor!(model, model.nodes, filter, monitor, :nodes)

    # available_vars = OrderedSet()
    # for node in monitor.nodes
    #     for dof in node.dofs
    #         push!(available_vars, dof.name)
    #         push!(available_vars, dof.natname)
    #     end
    # end
    # hint("NodeSumMonitor: ", replace(string(monitor.expr), " " => ""), " Available data: ", join(available_vars, ", "))
end


function update_monitor!(monitor::NodeSumMonitor, model)
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
            state[var] = evaluate(ex.args[2]; state...)
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
        monitor.vals[expr] = evaluate(expr; state...)
    end

    merge!(monitor.vals, extravals)
    monitor.vals[:stage] = model.env.stage
    monitor.vals[:T]     = model.env.T
    model.env.transient && (monitor.vals[:t]=model.env.t)
    push!(monitor.table, monitor.vals)

    # eval stop expressions
    for expr in monitor.stopexpr.args
        if evaluate(expr; state...)
            return failure("Stop condition at NodeSumMonitor ($expr)")
        end
    end

    return success()
end
