# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export @check

mutable struct AmaruException <: Exception
    message::String
end

Base.showerror(io::IO, e::AmaruException) = printstyled(io, "AmaruException: ", e.message, "\n", color=:red)

mutable struct PropertyException <: Exception
    message::String
end

Base.showerror(io::IO, e::PropertyException) = print(io, "PropertyException: ", e.message)

macro check(expr, exception=ArgumentError, args...)
    exprs = replace(string(expr), " "=>"")
    vars = getvars(expr)

    return quote
        if !$(esc(expr)) # Eval boolean expression
            if $exception==ArgumentError

                # Get function name
                st = stacktrace(backtrace())
                fname = :_
                for frame in st
                    if !startswith(string(frame.func), "_") && frame.func!=Symbol("macro expansion")
                        fname = frame.func
                        break
                    end
                end

                # Prepare message
                if length($vars)==1
                    prm = $(string(vars[1]))
                    val = $(esc(vars[1]))
                    msg = "Invalid value for $prm which must satisfy $($(exprs)). Got $prm = $val"
                else
                    prms = $(join(vars, ", ", " and "))
                    msg = "Arguments $prms must satisfy $($(exprs))"
                end

                fname != "" && (msg="$fname: $msg")
                throw($(exception)(msg, $(args)...))
            else
                throw($(exception)($(args)...))
            end
        end
    end
end

macro checkmissing(params, required, names)

    return quote
        # Get function name
        st = stacktrace(backtrace())
        fname = :_
        for frame in st
            if !startswith(string(frame.func), "_") && frame.func!=Symbol("macro expansion")
                fname = frame.func
                break
            end
        end

        missingkeys = setdiff($(esc(required)), keys($(esc(params))) )
        if length(missingkeys)>0
            msg = "Missing arguments: $(join(missingkeys, ", ")). Possible inputs are: $($(esc(names)))"
            msg = replace(msg, "=" => ":")
            error("$fname: $msg")
        end
    end
end