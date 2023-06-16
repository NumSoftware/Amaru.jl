# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export @check

mutable struct AmaruException <: Exception
    message::String
end

Base.showerror(io::IO, e::AmaruException) = print(io, "AmaruException: ", e.message)

macro check(expr, exception=ArgumentError, args...)
    exprs = replace(string(expr), " "=>"")
    vars = getvars(expr)

    return quote
        if !$(esc(expr)) # Eval boolean expression
            if $exception==ArgumentError

                # Get function name
                st = stacktrace(backtrace())
                fname = :_
                for frm in st
                    if !startswith(string(frm.func), "_") && frm.func!=Symbol("macro expansion")
                        fname = frm.func
                        break
                    end
                end

                # Prepare message
                if length($vars)==1
                    prm = $(string(vars[1]))
                    val = $(esc(vars[1]))
                    msg = "Invalid value for parameter $prm which must satisfy $($(exprs)). Got $prm = $val"
                else
                    prms = $(join(vars, ", ", " and "))
                    msg = "Parameters $prms must satisfy $($(exprs))"
                end

                fname != "" && (msg="$fname: $msg")
                throw($(exception)(msg, $(args)...))
            else
                throw($(exception)($(args)...))
            end
        end
    end
end