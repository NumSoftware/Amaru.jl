# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export Symbolic
export symbols, evaluate

mutable struct Symbolic
    sym::Symbol
    expr::Expr

    function Symbolic(symbol::Symbol)
        return new(symbol, :())
    end

    function Symbolic(expr::Expr)
        return new(:_, expr)
    end
end


# Gets current symbol or expression stored in a Symbolic object
function getexpr(symbolic::Symbolic)
    return symbolic.expr==:() ? symbolic.sym : symbolic.expr
end


getexpr(symbolic::Any) = symbolic


function Base.show(io::IO, symbolic::Symbolic)
    if symbolic.expr==:() 
        str = "Symbolic $(symbolic.sym)"
    else
        str = replace(string(symbolic.expr), r" ([*\^]) " => s"\1")
    end
    print(io, str)
end
 
for op in (:+, :-, :*, :/, :^, :>, :(>=), :<, :(<=), :(==), :(!=))
    @eval begin
        Base.$op(x::Symbolic, y::Symbolic) = Symbolic(Expr(:call, Symbol($op), getexpr(x), getexpr(y)))
        Base.$op(x::Symbolic, y) = Symbolic(Expr(:call, Symbol($op), getexpr(x), y))
        Base.$op(x, y::Symbolic) = Symbolic(Expr(:call, Symbol($op), x, getexpr(y)))
    end
end

for fun in (:abs, :sin, :cos, :tan, :log, :exp, :sqrt)
    @eval begin
        Base.$fun(x::Symbolic) = Symbolic(Expr(:call, Symbol($fun), getexpr(x)))
    end
end


for (fun, op) in zip((:and, :or), (:&&, :||))
    @eval begin
        $fun(x, y) = Symbolic(Expr($(QuoteNode(op)), getexpr(x), getexpr(y)))
        $fun(x, y, z) = $fun($fun(x,y), z)
        export $fun
    end
end


# macro symbols(syms...)
#     if syms[1] isa Expr && syms[1].head==:tuple 
#         syms = syms[1].args
#     end
#     expr = Expr(:block)
#     for symbolic in syms
#         push!(expr.args, :($(esc(symbolic)) = Symbolic($(QuoteNode(symbolic)))))
#     end

#     return expr
# end

macro define(syms...)

    expr = Expr(:block)
    for symbolic in syms
        if symbolic isa Symbol
            push!(expr.args, :($(esc(symbolic)) = Symbolic($(QuoteNode(symbolic)))))
        elseif symbolic isa Expr && symbolic.head == :(=) 
            lhs = symbolic.args[1]
            rhs = symbolic.args[2]
            !(lhs isa Symbol) && error("@define: Invalid definition $symbolic")
            push!(expr.args, :($(esc(lhs)) = $(esc(rhs))))
            push!(expr.args, :($(esc(lhs)).sym = $(QuoteNode(lhs))))
        else
            error("@define: Invalid definition $symbolic")
        end
    end
    push!(expr.args, nothing)
    return expr
end


@define x y z t
export @define, x, y, z, t

const arith_tol=1e-6

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
    :log10 => (a,) -> log10(a),
    :sqrt => (a,) -> sqrt(a),
    :max => (a,b) -> max(a,b),
    :min => (a,b) -> min(a,b),
    :length => (a,) -> length(a),
)


reduce!(arg; vars...) = arg


function reduce!(arg::Symbol; vars...)
    # should not be called on operators
    arg == :pi && return pi
    val = get(vars, arg, nothing)
    val === nothing && error("variable $arg not defined for this context")
    return val
end


function reduce!(expr::Expr; vars...)

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
        expr.args[i] = reduce!(expr.args[i]; vars...)
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

evaluate(expr::Real; vars...) = expr
evaluate(expr::Symbol; vars...) = reduce!(expr; vars...)
evaluate(expr::Expr; vars...) = reduce!(copy(expr); vars...)
evaluate(expr::Symbolic; vars...) = reduce!(copy(expr.expr); vars...)

@doc """
    evaluate(expr, vars...)

Returns the result of the arithmetic expression `expr` using values defined in `vars` if necessary.
""" evaluate


function getvars(symbolic::Symbol)
    if symbolic in ( :<, :<=, :>, >= )
        return []
    end
    return [ symbolic ]
end


function getvars(::Any)
    return []
end


function getvars(expr::Expr)
    symbols = Symbol[]
    expr.head==:. && error("getvars: dot operator not allowed")
    start = expr.head==:call ? 2 : 1
    for arg in expr.args[start:end]
        append!(symbols, getvars(arg))
    end
    return symbols
end


# replace a symbol
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


function fix_expr_maximum_minimum!(expr::Expr)
    if expr.head == :call && expr.args[1] in (:maximum, :minimum) && expr.args[2] isa Symbol
        symbol = expr.args[2]
        suffix = string(expr.args[1])[1:3]
        return Symbol("$(symbol)_$suffix")
    end

    start  = expr.head == :call ? 2 : 1
    finish = length(expr.args)

    for i in start:finish
        arg = expr.args[i]
        if isa(arg, Expr)
            expr.args[i] = fix_expr_maximum_minimum!(arg)
        end
    end
    return expr

end

function round_floats!(expr::Expr)
    start  = expr.head == :call ? 2 : 1
    finish = length(expr.args)

    for i in start:finish
        arg = expr.args[i]
        if arg isa AbstractFloat
            expr.args[i] = round(arg, digits=5)
        elseif arg isa Expr
            expr.args[i] = round_floats!(arg)
        end
    end
    return expr

end