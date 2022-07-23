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
    # @show expr
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
    # sym::Symbol
    expr::Expr

    function SymEx(symbol::Symbol)
        return new(Expr(:call, :identity, symbol))
    end

    function SymEx(expr::Expr)
        return new(expr)
    end

    # function SymEx(op::Symbol, args...)
    #     return new(op, collect(args))
    # end
end

function draw(a::SymEx) 
    return ifelse(a.expr.args[1] == :identity, a.expr.args[2], a.expr)
    # if a.expr.args[1] == :identity 
    #     return a.expr.args[2]
    # else
    #     return a.expr
    # end
end

draw(a) = a


for op in (:+, :-, :*, :/, :^, :>, :(>=), :<, :(<=), :(==), :(!=))
    @eval begin
        Base.$op(x::SymEx, y::SymEx) = SymEx(Expr(:call, Symbol($op), draw(x), draw(y)))
        Base.$op(x::SymEx, y) = SymEx(Expr(:call, Symbol($op), draw(x), y))
        Base.$op(x, y::SymEx) = SymEx(Expr(:call, Symbol($op), x, draw(y)))
    end
end

for fun in (:abs, :sin, :cos, :tan, :log, :exp)
    @eval begin
        Base.$fun(x::SymEx) = SymEx(Expr(:call, Symbol($fun), draw(x)))
    end
end



# function Base.show(io::IO, ex::SymEx)
#     nargs = length(ex.args)
#     if nargs==0
#         print(io, ex.root)
#     elseif nargs==1
#         op = ex.root
#         if op==:^
#             print(io, "(", ex.args[1], ",", op)
#         else
#             print(io, "(", ex.args[1], ex.root, ex.args[2], ")")
#         end
#     else
#         print(io, "(", ex.args[1], ex.root, ex.args[2], ")")
#     end
# end

# Base.:+(x::SymEx, y::SymEx) = SymEx(:+, x, y)
# Base.:+(x::SymEx, y) = SymEx(:+, x, y)
# Base.:+(x, y::SymEx) = SymEx(:+, x, y)
# Base.:-(x::SymEx, y::SymEx) = SymEx(:-, x, y)
# Base.:-(x::SymEx, y) = SymEx(:-, x, y)
# Base.:-(x, y::SymEx) = SymEx(:-, x, y)

# Base.:^(x::SymEx, y) = SymEx(:^, x, y)

# test(g) = @eval(Main, g(2))
# test(g) = eval(:(g(2)))
# test(g) = eval(Main, f(2)) #!


# Base.:+(x::Symbol, y::Symbol) = Expr(:call, :+, x, y)
# Base.:+(x::Expr, y::Expr) = Expr(:call, :+, x, y)
# Base.:+(x::Number, y::Symbol) = Expr(:call, :+, x, y)
# Base.:+(x::Symbol, y::Number) = Expr(:call, :+, x, y)
# Base.:+(x::Expr, y::Symbol) = Expr(:call, :+, x, y)
# Base.:+(x::Symbol, y::Expr) = Expr(:call, :+, x, y)

# for op in (:+, :-, :*, :/, :^, :>, :>=, :<, :<=)
#     @eval begin
#         Base.$op(x::Symbol, y::Symbol) = Expr(:call, Symbol($op), x, y)
#         Base.$op(x::Symbol, y) = Expr(:call, Symbol($op), x, y)
#         Base.$op(x, y::Symbol) = Expr(:call, Symbol($op), x, y)
#         Base.$op(x::Expr, y::Expr) = Expr(:call, Symbol($op), x, y)
#         Base.$op(x::Expr, y) = Expr(:call, Symbol($op), x, y)
#         Base.$op(x, y::Expr) = Expr(:call, Symbol($op), x, y)
#     end
# end

# Base.sin(x::Symbol) = Expr(:call, :sin, x)
# Base.cos(x::Symbol) = Expr(:call, :cos, x)
# Base.tan(x::Symbol) = Expr(:call, :tan, x)
# Base.exp(x::Symbol) = Expr(:call, :exp, x)
# Base.log(x::Symbol) = Expr(:call, :log, x)


# export Sym
# mutable struct Sym
#     sym::Symbol

#     function Sym(sym::Symbol)
#         return new(sym)
#     end
# end

# export SymEx
# mutable struct SymEx
#     expr::Expr

#     function SymEx(expr::Expr)
#         return new(expr)
#     end
# end

# draw(a::Sym) = a.sym
# draw(a::Expr) = a.expr
# draw(a) = a

# for op in (:+, :-, :*, :/, :^, :>, :(>=), :<, :(<=), :(==), :(!=))
#     @eval begin
#         # Base.$op(x::Sym, y::Sym) = SymEx(Expr(:call, Symbol($op), x.sym, y.sym))
#         # Base.$op(x::Sym, y) = SymEx(Expr(:call, Symbol($op), x.sym, draw(y)))
#         # Base.$op(x, y::Sym) = SymEx(Expr(:call, Symbol($op), draw(x), y.sym))

#         # Base.$op(x::SymEx, y::SymEx) = SymEx(Expr(:call, Symbol($op), x.expr, y.expr))
#         # Base.$op(x::SymEx, y) = SymEx(Expr(:call, Symbol($op), x.expr, draw(y)))
#         # Base.$op(x, y::SymEx) = SymEx(Expr(:call, Symbol($op), draw(x), y.expr))

#         Base.$op(x::SymEx, y::SymEx) = SymEx(Expr(:call, Symbol($op), x.expr, y.expr))
#         Base.$op(x::SymEx, y::Sym) = SymEx(Expr(:call, Symbol($op), x.expr, y.sym))
#         Base.$op(x::Sym, y::SymEx) = SymEx(Expr(:call, Symbol($op), x.sym, y.expr))
#         Base.$op(x::SymEx, y) = SymEx(Expr(:call, Symbol($op), x.expr, y))
#         Base.$op(x, y::SymEx) = SymEx(Expr(:call, Symbol($op), x, y.expr))
        
#         Base.$op(x::Sym, y::Sym) = SymEx(Expr(:call, Symbol($op), x.sym, y.sym))
#         Base.$op(x::Sym, y) = SymEx(Expr(:call, Symbol($op), x.sym, y))
#         Base.$op(x, y::Sym) = SymEx(Expr(:call, Symbol($op), x, y.sym))
#     end
# end





const x = SymEx(:x)
const y = SymEx(:y)
const z = SymEx(:z)
const t = SymEx(:t)
export x, y, z, t