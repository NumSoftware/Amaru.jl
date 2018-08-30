
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

function print_compact(io::IO, obj::Any)
    ty = typeof(obj) # update type

    if isbitstype(ty) || ty in [Symbol, Expr]
        print(io, obj)
    elseif ty<:AbstractString
        print(io, "\"", obj, "\"")
    elseif ty<:AbstractArray
        print_array_compact(obj)
    else
        print(io, ty)
    end
end


function print_array_compact(io::IO, array::AbstractArray)
    maxn = 8
    n = length(array)
    half = div(maxn,2)
    idx = n<=maxn ? [1:n;] : [1:half; n-half+1:n]
    print(io, "[")
    for i in idx
        print(io, array[i])
        if n>maxn && i==half
            print(io, ", … ")
        end
        if i!=n
            print(io, ", ")
        end
    end
    print(io, "]")
end


# Prints the values of an object fields
function print_field_values(io::IO, obj::Any)
    ctype = typeof(obj)
    print(io, ctype)
    for (field, ty) in zip(fieldnames(ctype), ctype.types )
        print(io, "\n  ", field, ": ")
        if !isdefined(obj, field)
            print(io, ty, "undef ")
            continue
        end

        item = getfield(obj, field)
        ty   = typeof(item) # update type


        if isbitstype(ty)
            print(io, item)
        elseif ty==String
            print(io, "\"", item, "\"")
        elseif :name in fieldnames(ty) 
            if isdefined(item, :name)
                print(io, ty, " name=", getfield(item, :name))
            else
                print(io, ty, " undef")
            end
        elseif :id in fieldnames(ty) 
            if isdefined(item, :id )
                print(io, ty, " id=", getfield(item, :id ))
            else
                print(io, ty, " undef")
            end
        #elseif ty <: AbstractVecOrMat{<:Real}
            #print(io, summary(item))
        elseif ty<:AbstractArray
            s = summary(item)
            if length(s)>40
                s = s[ 1: findfirst("{",s).start-1 ]
            end
            print(io, s)
        elseif ty<:Function
            print(io, "Function ", item)
        else
            print(io, summary(item))
        end
    end
    return nothing
end

# Print arrays of objects
function print_array_values(io::IO, array::AbstractArray)
    print(io, summary(array))
    n = length(array)
    n==0 && return
    print(io, ":")
    maxn = 8
    half = div(maxn,2)
    idx = n<=maxn ? [1:n;] : [1:half; n-half+1:n]
    for i in idx
        print(io, "\n  ", i, ": ", replace(string(array[i]), "\n" => "\n  "))
        if n>maxn && i==half
            print(io, "\n    ⋮")
        end
    end
    return
end
