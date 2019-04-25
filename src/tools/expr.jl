# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

const arith_tol=1e-6

const op_dict = Dict{Symbol,Function}(
    :+ => +, :- => -, :* => *, :/ => /, :^ => ^, 
    :div => div,
    :(>)  => (a,b) -> a>b+arith_tol,
    :(<)  => (a,b) -> a<b-arith_tol,
    :(>=) => (a,b) -> a>b-arith_tol,
    :(<=) => (a,b) -> a<b+arith_tol,
    :(==) => (a,b) -> abs(a-b)<arith_tol)


export eval_arith_expr

reduce_arith_expr!(arg; vars...) = arg

function reduce_arith_expr!(arg::Symbol; vars...)
    # should not be called on operators
    val = get(vars, arg, nothing)
    val == nothing && error("variable $arg not defined for this context")
    return val
end


function reduce_arith_expr!(expr::Expr; vars...)

    # geto operator arguments
    #op_args = expr.head==:call ? expr.args[2:end] : expr.args
    arg_idxs = expr.head==:call ? UnitRange(2,length(expr.args)) : UnitRange(1,length(expr.args))

    # recursively reduce on operator arguments and replace variables
    #for i=1:length(op_args)
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
    end

    return expr
end


eval_arith_expr(expr::Real; vars...) = expr
eval_arith_expr(expr::Expr; vars...) = reduce_arith_expr!(copy(expr); vars...)
eval_arith_expr(expr::Symbol; vars...) = reduce_arith_expr!(copy(expr); vars...)


# Fixes comparisons expressions using a tolerance
function fix_comparison_scalar(expr::Expr)
    mexpr = copy(expr) # expression to be modified
    tol = 1e-6

    fix_comp = function(expr::Expr)
        symb = expr.args[1]
        a = expr.args[2]
        b = expr.args[3]
        if typeof(a)==String || typeof(b)==String
            return expr
        end
        if symb == :(==)
            return :(abs($a-$b) < $tol)
        end
        if symb == :(>=)
            return :($a > $b - $tol)
        end
        if symb == :(>)
            return :($a > $b + $tol)
        end
        if symb == :(<=)
            return :($a < $b + $tol)
        end
        if symb == :(<)
            return :($a < $b - $tol)
        end
        return expr
    end

    if mexpr.head == :call
        return fix_comp(mexpr)
    else
        for (i,arg) in enumerate(mexpr.args)
            if typeof(arg)!=Expr; continue end
            if arg.head == :call
                mexpr.args[i] = fix_comp(arg)
            else
                mexpr.args[i] = fix_comparison_scalar(arg)
            end
        end
    end
    return mexpr
end


# Define a function based on an expression
export @fun
export @cond

macro fun(expr)
    return quote
        (t, x, y, z) -> $expr
    end
end

macro cond(expr)
    expr = fix_comparison_scalar(expr)
    return quote
        (t, x, y, z) -> $expr
    end
end


# Fixes comparisons in array expressions using a tolerance
function fix_comparison_arrays(expr::Expr)
    mexpr = copy(expr) # expression to be modified
    tol = 1e-6

    fix_comp = function(expr::Expr)
        symb = expr.args[1]
        a = expr.args[2]
        b = expr.args[3]
        if typeof(a)==String || typeof(b)==String
            return expr
        end
        if symb == :(==)
            return :(maximum(abs.($a.-$b)) < $tol)
        end
        if symb == :(>=)
            return :(minimum($a) > maximum($b) - $tol)
        end
        if symb == :(>)
            return :(minimum($a) > maximum($b) + $tol)
        end
        if symb == :(<=)
            return :(maximum($a) < minimum($b) + $tol)
        end
        if symb == :(<)
            return :(maximum($a) < minimum($b) - $tol)
        end
        return expr
    end

    if mexpr.head == :call && length(mexpr.args)==3 # check comparison
        return fix_comp(mexpr)
    else
        for (i,arg) in enumerate(mexpr.args)
            if typeof(arg)!=Expr; continue end
            if arg.head == :call
                mexpr.args[i] = fix_comp(arg)
            else
                mexpr.args[i] = fix_comparison_arrays(arg)
            end
        end
    end
    return mexpr
end
