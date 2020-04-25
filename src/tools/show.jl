
# show the content of an array variable or expression including calling file
macro showm(M)
    filename = @__FILE__
    return quote
        printstyled($(string(M)), " = \n", color=:light_magenta)
        #Base.showarray(STDOUT, $(esc(M)), false)
        Base.print_matrix(stdout, $(esc(M)))
        #filename != nothing && printstyled(" : ", basename($filename), color=:blue)
        println()
    end
end

macro display(x)
    return quote
        printstyled($(string(x)), " = \033[K\n", color=:light_magenta)
        display($(esc(x)))
    end
end


export custom_dump


custom_dump(x; maxdepth=2) = custom_dump(stdout, x, maxdepth, "")


function custom_dump(io::IO, x, maxdepth::Int, indent)
    T = typeof(x)

    nf = nfields(x)

    if nf==0 
        print(io, repr(x))
        return
    end

    print(io, T)

    if maxdepth==0
        if :name in fieldnames(T) && isdefined(x, :name )
            print(io, "  name=", repr(getfield(x, :name)))
        end
        if :id in fieldnames(T) && isdefined(x, :id )
            print(io, "  id=", getfield(x, :id))
        end
        if :tag in fieldnames(T) && isdefined(x, :tag )
            tag = getfield(x, :tag)
            tag=="" || print(io, "  tag=", repr(tag))
        end
        return
    end

    maxdepth==0 && return

    for field in 1:nf
        fname = string(fieldname(T, field))
        fname[1]=='_' && continue
        println(io)
        print(io, indent, "  ", fname, ": ")
        if isdefined(x,field)
            custom_dump(io, getfield(x, field), maxdepth-1, indent*"  ")
        else
            print(io, "#undef")
        end
    end

    return

end


custom_dump(io::IO, x::Function, maxdepth::Int, indent) = print(io, "Function ", x)


function custom_dump(io::IO, x::AbstractDict, maxdepth::Int, indent)
    s = split(summary(x), ".")[end]
    #print(io, summary(x))
    print(io, s)
    maxdepth==0 && return

    for (k,v) in x
        println(io)
        print(io, indent*"  ", repr(k), " => ")
        custom_dump(io, v, 0, indent*"  ")
    end
end


# Print arrays of objects
function custom_dump(io::IO, x::AbstractArray, maxdepth::Int, indent="")
    n  = length(x)

    if eltype(x) <: Number
        if n<20 
            s = repr(x)
            if length(s)<70
                show(io, x)
                return
            end
        end
        print(io, summary(x))
        return
    end

    print(io, summary(x))

    maxdepth==0 && return

    print(io, ":")
    maxn = 8
    half = div(maxn,2)
    idx = n<=maxn ? [1:n;] : [1:half; n-half+1:n]
    for i in idx
        println(io)
        print(io, indent*"  ", i, ": ")
        ety = eltype(x)
        if nfields(x[i])==0
            show(io, x[i])
        else
            custom_dump(io, x[i], maxdepth-1, indent*"  ")
        end
        if n>maxn && i==half
            println(io)
            print(io, indent*"  ", "⋮")
        end
    end
    return
end
