# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

const arith_tol=1e-6

const op_dict = Dict{Symbol,Function}(
    :+ => +,
    :- => -,
    :* => *,
    :/ => /,
    :^ => ^,
    :div => div,
    :(>)  => (a,b) -> a>b+arith_tol,
    :(<)  => (a,b) -> a<b-arith_tol,
    :(>=) => (a,b) -> a>b-arith_tol,
    :(<=) => (a,b) -> a<b+arith_tol,
    :(==) => (a,b) -> abs(a-b)<arith_tol,
    :(!=) => (a,b) -> abs(a-b)>=arith_tol,
    :abs => (a,) -> abs(a),
    :sin => (a,) -> sin(a),
    :cos => (a,) -> cos(a),
    :tan => (a,) -> tan(a),
    :exp => (a,) -> exp(a),
    :log => (a,) -> log(a),
   )


export eval_arith_expr

reduce_arith_expr!(arg; vars...) = arg

function reduce_arith_expr!(arg::Symbol; vars...)
    # should not be called on operators
    arg == :pi && return pi
    val = get(vars, arg, nothing)
    val === nothing && error("variable $arg not defined for this context")
    return val
end


function reduce_arith_expr!(expr::Expr; vars...)

    # get operation arguments
    if expr.head == :call
        arg_idxs = UnitRange(2,length(expr.args))
    elseif expr.head == :comparison
        arg_idxs = StepRange(1,2,length(expr.args))
    else # &&, ||, ?:, comparison
        arg_idxs = UnitRange(1,length(expr.args))
    end

    # recursively reduce on arguments and replace variables
    for i in arg_idxs
        expr.args[i] = reduce_arith_expr!(expr.args[i]; vars...)
    end

    # evaluate operations
    if expr.head == :call
        op = get(op_dict, expr.args[1], nothing)
        op == nothing && error("operation $(expr.args[1]) not allowed in this context")
        return op(expr.args[2:end]...)
    elseif expr.head == :&&
        return expr.args[1] && expr.args[2]
    elseif expr.head == :||
        return expr.args[1] || expr.args[2]
    elseif expr.head == :comparison
        for i=2:2:length(expr.args)
            op = op_dict[expr.args[i]]
            op(expr.args[i-1], expr.args[i+1]) || return false
        end
        return true
    elseif expr.head == :if
        return expr.args[1] ? expr.args[2] : expr.args[3]
    end

    return expr
end

eval_arith_expr(expr::Real; vars...) = expr
eval_arith_expr(expr::Expr; vars...) = reduce_arith_expr!(copy(expr); vars...)
eval_arith_expr(expr::Symbol; vars...) = reduce_arith_expr!(expr; vars...)

@doc """
    eval_arith_expr(expr, vars...)

Returns the result of the arithmetic expression `expr` using values defined in `vars` if necessary.
""" eval_arith_expr


function get_vars(sym::Symbol)
    return [ sym ]
end

function get_vars(sym::Any)
    return []
end

function get_vars(expr::Expr)
    symbols = Symbol[]
    start = expr.head==:call ? 2 : 1
    for arg in expr.args[start:end]
        append!(symbols, get_vars(arg))
    end
    return symbols
end

export @check
macro check(expr)
    exprs = replace(string(expr), " "=>"")
    vars = get_vars(expr)

    return quote
        if !$(esc(expr)) # Eval boolean expression
            # Get function name
            st = stacktrace(backtrace())
            fname = ""
            for frm in st
                if frm.func != :backtrace && frm.func!= Symbol("macro expansion")
                    fname = frm.func
                    break
                end
            end

            # Prepare message
            if length($vars)==1
                prm = $(string(vars[1]))
                val = $(esc(vars[1]))
                msg = "Invalid value for parameter $prm ($val) which must satisfy $($(exprs))"
            else
                prms = $(join(vars, ", ", " and "))
                msg = "Parameters $prms must satisfy $($(exprs))"
            end

            fname != "" && (msg="$fname: $msg")
            error(msg)
        end
    end
end


export SymEx
mutable struct SymEx
    root::Symbol
    args::Array{Any,1}

    function SymEx(symbol::Symbol)
        return new(symbol, [])
    end

    function SymEx(op::Symbol, args...)
        return new(op, collect(args))
    end
end

function Base.show(io::IO, ex::SymEx)
    nargs = length(ex.args)
    if nargs==0
        print(io, ex.root)
    elseif nargs==1
        op = ex.root
        if op==:^
            print(io, "(", ex.args[1], ",", op)
        else
            print(io, "(", ex.args[1], ex.root, ex.args[2], ")")
        end
    else
        print(io, "(", ex.args[1], ex.root, ex.args[2], ")")
    end
end

Base.:+(x::SymEx, y::SymEx) = SymEx(:+, x, y)
Base.:+(x::SymEx, y) = SymEx(:+, x, y)
Base.:+(x, y::SymEx) = SymEx(:+, x, y)
Base.:-(x::SymEx, y::SymEx) = SymEx(:-, x, y)
Base.:-(x::SymEx, y) = SymEx(:-, x, y)
Base.:-(x, y::SymEx) = SymEx(:-, x, y)

Base.:^(x::SymEx, y) = SymEx(:^, x, y)

test(g) = @eval(Main, g(2))
# test(g) = eval(:(g(2)))
# test(g) = eval(Main, f(2))