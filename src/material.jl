# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

import Base.copy!

# Abstract type for material


"""
`MatParams`

Abstract type for objects used to store material parameters and to define
the behaviour of elements.
"""
abstract type MatParams end
Base.:(<<)(a::Tuple, b::MatParams) = return (a..., b)
Base.:(<<)(a, b::MatParams) = return (a, b)
Base.:(=>)(a::Tuple, b::MatParams) = return (a.first, a.second, b)
Base.:(=>)(a, b::MatParams) = return (a, b)



# function copy!(target::MatParams, source::MatParams)
#     @show("Trial function: This sould not be used")
#     # source and target must be of the same type
#     T = typeof(source)
#     for fld in fieldnames(source)
#         if fieldtype(T, fld) <: Array
#             if isdefined(target,fld)
#                 getfield(target, fld) .= getfield(source, fld)
#             else
#                 setfield!(target, fld, copy(getfield(source, fld)))
#             end
#         else
#             setfield!(target, fld, getfield(source, fld))
#         end
#     end
# end

"""
`output_keys(matparams)`

Returns a list of keys from an specified material
for output at element level or nodal level.
"""
function output_keys(matparams::MatParams)
    return Symbol[]
end


function paramsdict(matparams::MatParams)
    return OrderedDict( string(field)=> getfield(matparams, field) for field in fieldnames(typeof(matparams)) )
end


function databook(mats::Array{<:Pair,1})
    db = DataBook()
    for (label, matparams) in mats
        # tab = DataTable(;paramsdict(matparams)...)
        dt = DataTable()
        push!(dt, paramsdict(matparams))
        setfield!(dt, :name, string(label))
        push!(db, dt)
    end
    return db
end