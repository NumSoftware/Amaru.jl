# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export Symbolic
mutable struct Symbolic
    expr::Expr

    function Symbolic(symbol::Symbol)
        return new(Expr(:call, :identity, symbol))
    end

    function Symbolic(expr::Expr)
        return new(expr)
    end
end

function getexpr(a::Symbolic)
    return a.expr.args[1]==:identity ? a.expr.args[2] : a.expr
end

getexpr(a) = a

function Base.show(io::IO, ex::Symbolic)
    expr = getexpr(ex)
    str = replace(string(expr), r" ([*\^]) " => s"\1")
    print(io, str)
end
 
for op in (:+, :-, :*, :/, :^, :>, :(>=), :<, :(<=), :(==), :(!=))
    @eval begin
        Base.$op(x::Symbolic, y::Symbolic) = Symbolic(Expr(:call, Symbol($op), getexpr(x), getexpr(y)))
        Base.$op(x::Symbolic, y) = Symbolic(Expr(:call, Symbol($op), getexpr(x), y))
        Base.$op(x, y::Symbolic) = Symbolic(Expr(:call, Symbol($op), x, getexpr(y)))
    end
end

for fun in (:abs, :sin, :cos, :tan, :log, :exp)
    @eval begin
        Base.$fun(x::Symbolic) = Symbolic(Expr(:call, Symbol($fun), getexpr(x)))
    end
end

for fun in (:and, :or)
    @eval begin
        $fun(x::Symbolic, y::Symbolic) = Symbolic(Expr(:&&,  getexpr(x), getexpr(y)))
        export $fun
    end
end

function Base.convert(::Type{Expr}, ex::Symbolic)
    return ex.expr
end

macro vars(symbols...)
    expr = Expr(:block)
    for sym in symbols
        push!(expr.args, :($(esc(sym)) = Symbolic($(QuoteNode(sym)))))
    end

    return expr
end
export @vars


const arith_tol=1e-5

const op_dict = Dict{Symbol,Function}(
    :+ => +,
    :- => -,
    :* => *,
    :/ => /,
    :^ => ^,
    :and => (a,b) -> a && b,
    :or => (a,b) -> a || b,
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
        op === nothing && error("operation $(expr.args[1]) not allowed in this context")
        return op(expr.args[2:end]...)
    elseif expr.head == :&&
        return expr.args[1] && expr.args[2]
    elseif expr.head == :||
        return expr.args[1] || expr.args[2]
    elseif expr.head == :comparison
        for i in 2:2:length(expr.args)
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
eval_arith_expr(expr::Symbol; vars...) = reduce_arith_expr!(expr; vars...)
eval_arith_expr(expr::Expr; vars...) = reduce_arith_expr!(copy(expr); vars...)
eval_arith_expr(expr::Symbolic; vars...) = reduce_arith_expr!(copy(expr.expr); vars...)

@doc """
    eval_arith_expr(expr, vars...)

Returns the result of the arithmetic expression `expr` using values defined in `vars` if necessary.
""" eval_arith_expr


function getvars(sym::Symbol)
    if sym in ( :<, :<=, :>, >= )
        return []
    end
    return [ sym ]
end

function getvars(::Any)
    return []
end

function getvars(expr::Expr)
    symbols = Symbol[]
    start = expr.head==:call ? 2 : 1
    for arg in expr.args[start:end]
        append!(symbols, getvars(arg))
    end
    return symbols
end

function Base.replace(expr::Expr, pair::Pair)
    src    = pair.first
    tgt    = pair.second
    expr   = copy(expr)
    start  = expr.head == :call ? 2 : 1
    finish = length(expr.args)

    for i in start:finish
        arg = expr.args[i]
        if isa(arg, Symbol) && arg==src
            expr.args[i] = tgt
        elseif isa(arg, Expr)
            expr.args[i] = replace(arg, pair)
        end
    end
    return expr
end


# const x = Symbolic(:x)
# const y = Symbolic(:y)
# const z = Symbolic(:z)
# const t = Symbolic(:t)
# export x, y, z, t