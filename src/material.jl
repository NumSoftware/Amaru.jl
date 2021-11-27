# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

import Base.copy!

# Abstract type for material
# ==========================

"""
`Material`

Abstract type for objects used to store material parameters and to define
the behaviour of elements.
"""
abstract type Material end

function copy!(target::Material, source::Material)
    @show("Trial function: This sould not be used")
    # source and target must be of the same type
    T = typeof(source)
    for fld in fieldnames(source)
        if fieldtype(T, fld) <: Array
            if isdefined(target,fld)
                getfield(target, fld) .= getfield(source, fld)
            else
                setfield!(target, fld, copy(getfield(source, fld)))
            end
        else
            setfield!(target, fld, getfield(source, fld))
        end
    end
end

"""
`output_keys(mat)`

Returns a list of keys from an specified material
for output at element level or nodal level.
"""
function output_keys(mat::Material)
    return Symbol[]
end
