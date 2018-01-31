# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

struct Functor
    rhs::Union{Real, Symbol, Expr}
    fun::Function
    cval::Real # constant value
    function Functor(argnames::Union{Symbol,Expr}, rhs::Union{Real, Symbol, Expr})
        try
            fun = eval( :($argnames -> $rhs) )
            typeof(rhs)<:Real && return new(rhs, fun, rhs)
            return new(rhs, fun)
        catch ex
            println("Functor: cannot create Functor with parameters :($argnames) and :($rhs)")
            throw(ex)
        end
    end
end

function (functor::Functor)(args...)
    isdefined(functor, :cval) && return functor.cval
    return Base.invokelatest(functor.fun, args...)
end
