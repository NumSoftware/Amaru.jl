# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

import Base.copy!

# Abstract type for material


"""
`Material`

Abstract type for objects used to store material parameters and to define
the behaviour of elements.
"""
abstract type Material end
# Base.:(<<)(a::Tuple, b::Material) = return (a..., b)
# Base.:(<<)(a, b::Material) = return (a, b)
# Base.:(=>)(a::Tuple, b::Material) = return (a.first, a.second, b)
# Base.:(=>)(a, b::Material) = return (a, b)

# Base.:(<<)(a::Tuple, b::Type{<:Material}) = return (a..., b)
# Base.:(<<)(a, b::Type{<:Material}) = return (a, b)
# Base.:(<<)(a::Type{<:Material}, b) = return (a, b)


# function copy!(target::Material, source::Material)
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
`output_keys(mat)`

Returns a list of keys from an specified material
for output at element level or nodal level.
"""
function output_keys(mat::Material)
    return Symbol[]
end


function paramsdict(mat::Material)
    return OrderedDict( string(field)=> getfield(mat, field) for field in fieldnames(typeof(mat)) )
end


function databook(mats::Array{<:Pair,1})
    db = DataBook()
    for (label, mat) in mats
        # tab = DataTable(;paramsdict(mat)...)
        dt = DataTable()
        push!(dt, paramsdict(mat))
        setfield!(dt, :name, string(label))
        push!(db, dt)
    end
    return db
end