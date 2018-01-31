# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

import Base.show
export show

function print_short(io::IO, obj::Any)
    ty = typeof(obj) # update type
    if isbits(ty) || ty in [Symbol, Expr]
        print(io, obj)
    elseif ty<:AbstractString
        print(io, "\"", obj, "\"")
    elseif ty<:AbstractArray
        print_short_array(obj)
    #elseif ty<:Dict
        #print_short_dict(obj)
    else
        print(io, ty)
    end
end


function print_short_array(io::IO, array::AbstractArray)
    const maxn = 8
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
function print_datatype_fields(io::IO, obj::Any)
    ctype = typeof(obj)
    print(io, ctype, " with:")
    for (field, ty) in zip(fieldnames(ctype), ctype.types )
        print(io, "\n  ", field, ": ")
        if !isdefined(obj, field)
            print(io, ty, "#undef ")
            continue
        end

        item = getfield(obj, field)
        ty   = typeof(item) # update type
        if :name in fieldnames(ty) 
            if isdefined(item, :name)
                print(io, ty, ", name field ", getfield(item, :name))
            end
        elseif :id in fieldnames(ty) 
            if isdefined(item, :id )
                print(io, ty, ", id field ", getfield(item, :id ))
            end
        elseif ty<:AbstractArray && ndims(ty)==1
            array = getfield(obj, field)
            n = length(array)
            if eltype(ty)<:Real && length(size(array))==1
                print_short_array(io, array)
            else
                if n==0
                    print(io, ty, " [ ]")
                else
                    print(io, ty, ", length ", length(array))
                end
            end

            if n>0
                maxn = 10
                # add items IDs if available
                has_id = all( :id in fieldnames(item) for item in array[1:min(maxn,n)] )
                if has_id
                    ids = [ item.id for item in array ]
                    print(io, ", id fields ")
                    print_short_array(io, ids)
                end
                # add items names if available
                has_name = all( :name in fieldnames(item) for item in array[1:min(maxn,n)] )
                if has_name
                    names = [ item.name for item in array ]
                    print(io, ", name fields ")
                    print_short_array(io, names)
                end
            end
        elseif ty<:AbstractArray
            array = getfield(obj, field)
            str   = replace(string(size(array)), ", ", "×")[2:end-1]
            print(io, ty, ", size ", str)
        else
            print_short(io, item)
        end
    end
    return nothing
end

# Print arrays of objects
function print_datatype_array(io::IO, array::AbstractArray)
    print(io, length(array), "-element ",typeof(array), ":")
    n = length(array)
    maxn = 8
    half = div(maxn,2)
    idx = n<=maxn ? [1:n;] : [1:half; n-half+1:n]
    for i in idx
        str = "\n  "*replace(string(array[i]), "\n","\n  ")
        print(io, str)
        if n>maxn && i==half
            print(io, "\n    ⋮")
        end
    end
    return nothing
end


# Generates show functions for data types
macro show_function(datatype)
    return quote
        function $(esc(:show))(io::IO, obj::$datatype)
            print_datatype_fields(io, obj)
        end
    end
end

# Reuses show function to display array of objects
macro show_array_function(datatype)
    return quote
        function $(esc(:show))(io::IO, array::Array{<:$(datatype),1})
            print_datatype_array(io, array)
        end
    end
end

