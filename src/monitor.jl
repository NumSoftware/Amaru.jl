# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# Log structs and functions
# =========================

abstract type AbstractMonitor end


# Ip Monitor
# ==========

mutable struct IpMonitor<:AbstractMonitor
    expr     ::Union{Symbol,Expr}
    val      ::Number
    filename ::String
    filter   ::Union{Int,Symbol,String,Expr}
    table    ::DataTable
    ip       ::Ip

    function IpMonitor(expr::Union{Symbol,Expr}, filename::String="")
        return new(expr, NaN, filename, :(), DataTable())
    end
end

# function output(mon::IpMonitor) 
#     val = round(mon.val, sigdigits=5)
#     return "ip $(mon.ip.id)  $(mon.expr) $val"
# end


function setup_monitor!(domain, filter, monitor::IpMonitor)
    monitor.filter = filter
    if filter isa Integer
        monitor.ip = domain.elems[:ips][filter]
        return
    end
    ips = domain.elems[:ips][filter]
    n = length(ips)
    n == 0 && warn("setup_monitor!: No ips found for filter expression: $(monitor.filter)")
    n >  1 && notify("setup_monitor!: More than one ip match filter expression: $(monitor.filter)")
    n >= 1 && (monitor.ip = ips[1])
    return nothing
end


function update_monitor!(monitor::IpMonitor, domain)
    isdefined(monitor, :ip) || return

    state       = ip_vals(monitor.ip)
    monitor.val = eval_arith_expr(monitor.expr; state...)

    # data = OrderedDict(:val=>monitor.val)
    # data = OrderedDict(:stage=> domain.env.cstage, :T=>domain.env.T)
    # domain.env.transient && (vals[:t] = domain.env.t)

    monitor.vals[:stage] = domain.env.cstage
    monitor.vals[:T]     = domain.env.T
    domain.env.transient && (monitor.vals[:t]=dom.env.t)
    push!(monitor.table, data)

    if monitor.filename!="" 
        filename = joinpath(domain.env.outdir, monitor.filename)
        save(monitor.table, filename, printlog=true)
    end
end

# Node Monitor
# ============

mutable struct NodeMonitor<:AbstractMonitor
    expr     ::Union{Symbol,Expr}
    vals     ::OrderedDict # stores current values
    filename ::String
    filter   ::Union{Int,Symbol,String,Expr}
    table    ::DataTable
    node     ::Node

    function NodeMonitor(expr::Union{Symbol,Expr}, filename::String="")
        if expr isa Symbol || expr.head == :call
            expr = :($expr,)
        end

        return new(expr, OrderedDict(), filename, :(), DataTable())
    end
end

# function output(monitor::NodeMonitor) 
#     str = "node $(replace(string(monitor.filter), " " => ""))"
#     for (k,v) in monitor.vals
#         k in (:stage, :T) && continue
#         v = v isa Real ? round(v, sigdigits=5) : v
#         str *= "  $(replace(string(k), " " => "")):$v" 
#     end

#     return str
# end

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


function setup_monitor!(domain, filter, monitor::NodeMonitor)
    monitor.filter = filter
    if filter isa Integer
        monitor.node = domain.elems[:nodes][filter]
        return
    end
    nodes = domain.elems[:nodes][filter]
    n = length(nodes)
    n == 0 && warn("setup_monitor!: No nodes found for filter expression: $(monitor.filter)")
    n >  1 && notify("setup_monitor!: More than one node match filter expression: $(monitor.filter)")
    n >= 1 && (monitor.node = nodes[1])
    return nothing
end


function update_monitor!(monitor::NodeMonitor, domain)
    isdefined(monitor, :node) || return

    for expr in monitor.expr.args
        state = node_vals(monitor.node) 
        monitor.vals[expr] = eval_arith_expr(expr; state...)
    end

    monitor.vals[:stage] = domain.env.cstage
    monitor.vals[:T]     = domain.env.T
    domain.env.transient && (monitor.vals[:t]=dom.env.t)

    push!(monitor.table, monitor.vals)

    if monitor.filename!=""
        filename = joinpath(domain.env.outdir, monitor.filename)
        save(monitor.table, filename, printlog=true)
    end
end



# Monitor for a group of ips
# ==========================

mutable struct IpGroupMonitor<:AbstractMonitor
    expr    ::Expr
    vals    ::OrderedDict
    filename::String
    filter  ::Union{Symbol,String,Expr}
    table   ::DataTable
    ips     ::Array{Ip,1}

    function IpGroupMonitor(expr::Expr, filename::String="")
        if expr.head == :call
            expr = :($expr,)
        end

        return new(expr, OrderedDict(), filename, :(), DataTable(), [])
    end
end


# function output(monitor::IpGroupMonitor)
#     str = "ip group $(replace(string(monitor.filter), " " => ""))"
#     for (k,v) in monitor.vals
#         k in (:stage, :T) && continue
#         str *= "  $(replace(string(k), " " => "")):$v" 
#     end

#     return str
# end


function setup_monitor!(domain, filter, monitor::IpGroupMonitor)
    monitor.filter = filter
    monitor.ips = domain.elems[:ips][monitor.filter]
    length(monitor.ips)==0 && warn("setup_monitor!: No ips found for filter expression: ", monitor.filter)
end


function update_monitor!(monitor::IpGroupMonitor, domain)
    length(monitor.ips) == 0 && return

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

    monitor.vals[:stage] = domain.env.cstage
    monitor.vals[:T]     = domain.env.T
    domain.env.transient && (monitor.vals[:t]=dom.env.t)

    push!(monitor.table, monitor.vals)

    if monitor.filename!="" 
        filename = joinpath(domain.env.outdir, monitor.filename)
        save(monitor.table, filename, printlog=true)
    end
end