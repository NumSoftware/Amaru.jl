
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
        printstyled($(string(x)), " = \e[K\n", color=:light_magenta)
        display($(esc(x)))
    end
end

macro s(exs...)
    blk = Expr(:block)
    sf = basename(string(__source__.file))
    sl = string(__source__.line)
    for ex in exs
        push!(blk.args, 
        :(
            printstyled($(sprint(Base.show_unquoted,ex)*" = "),
            repr(begin value=$(esc(ex)) end), color=13);
            printstyled("  ", $sf, ":", $sl, "\n", color=:light_black)
        ))
    end
    isempty(exs) || push!(blk.args, :value)
    return blk
end


export custom_dump


custom_dump(x; maxdepth=2) = custom_dump(stdout, x, maxdepth, "")

not_unrolling_typenames_in_custom_dump = [ "ShapeType" ]

function custom_dump(io::IO, x, maxdepth::Int, indent::String)
    T = typeof(x)

    nf = nfields(x)

    # Prints objects with no fields
    if nf==0 
        print(io, repr(x))
        return
    end

    print(io, T)

    if maxdepth==0
        if :name in fieldnames(T) && isdefined(x, :name)
            print(io, "  name=", repr(getfield(x, :name)))
        end
        if :id in fieldnames(T) && isdefined(x, :id )
            print(io, "  id=", getfield(x, :id))
        end
        if :tag in fieldnames(T) && isdefined(x, :tag )
            tag = getfield(x, :tag)
            tag=="" || print(io, "  tag=", repr(tag))
        end
    end

    maxdepth==0 && return

    for field in 1:nf

        fname = string(fieldname(T, field))
        fname[1]=='_' && continue
        println(io)


        print(io, indent, "  ", fname, ": ")


        if isdefined(x, field)
            fieldvalue = getfield(x, field)

            if string(typeof(fieldvalue)) in not_unrolling_typenames_in_custom_dump
                custom_dump(io, fieldvalue, 0, indent*"  ")
                continue
            end

            custom_dump(io, fieldvalue, maxdepth-1, indent*"  ")
        else
            print(io, "#undef")
        end
    end

    return

end


custom_dump(io::IO, x::Function, maxdepth::Int, indent::String) = print(io, "Function ", x)


function custom_dump(io::IO, x::AbstractDict, maxdepth::Int, indent::String)
    if valtype(x)<:Real || valtype(x)<:AbstractString
        n = length(x)
        if n<8 
            s = repr(x)
            if length(s)<70
                show(io, x)
                return
            end
        end
    end

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
function custom_dump(io::IO, x::AbstractArray, maxdepth::Int, indent::String="")
    n = length(x)
    if n==0
        print(io, summary(x))
        return
    end

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
