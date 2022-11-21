
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


export _show

_show(x; maxdepth=2) = _show(stdout, x, maxdepth, "", "    ")

const not_unrolling_typenames_in_show = [ "CellShape", "ModelEnv"  ]

function _show(io::IO, x, maxdepth::Int, indent::String="", tab::String="    ")
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


        print(io, indent*tab, fname, ": ")


        if isdefined(x, field)
            fieldvalue = getfield(x, field)

            if string(typeof(fieldvalue)) in not_unrolling_typenames_in_show
                _show(io, fieldvalue, 0, indent*tab, tab)
                continue
            end

            _show(io, fieldvalue, maxdepth-1, indent*tab)
        else
            print(io, "#undef")
        end
    end

    return

end


_show(io::IO, x::Function, maxdepth::Int, indent::String, tab::String) = print(io, "Function ", x)
_show(io::IO, x::Expr, maxdepth::Int, indent::String, tab::String) = print(io, "Expr ", x)


function _show(io::IO, x::AbstractDict, maxdepth::Int, indent::String="", tab::String="    ")
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
        print(io, indent*tab, repr(k), " => ")
        _show(io, v, 0, indent*tab)
    end
end


# Print arrays of objects
function _show(io::IO, x::Union{AbstractArray, AbstractSet}, maxdepth::Int, indent::String="", tab::String="    ")
    summ = summary(x)
    n = length(x)
    if n==0
        print(io, summ)
        return
    end

    isset = x isa AbstractSet
    isset && (x = collect(x))

    if eltype(x) <: Number
        if n<20 
            s = repr(x)
            if length(s)<70
                show(io, x)
                return
            end
        end
        print(io, summ)
        return
    end

    print(io, summ)

    maxdepth==0 && return

    print(io, ":")
    maxn = 20
    half = div(maxn,2)
    idx = n<=maxn ? [1:n;] : [1:half; n-half+1:n]
    for i in idx
        println(io)
        print(io, indent*tab, i, ": ")
        ety = eltype(x)
        if nfields(x[i])==0
            show(io, x[i])
        else
            _show(io, x[i], maxdepth-1, indent*tab)
        end
        if n>maxn && i==half
            println(io)
            print(io, indent*tab, "⋮")
        end
    end
end
