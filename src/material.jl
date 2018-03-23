# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

import Base.copy!

# Abstract type for material
# ==========================

"""
`Material`

Abstract type for objects used to store material parameters and to define
the bejaviour of elements.
"""
abstract type Material end

function copy!(target::Material, source::Material)
    @show("Trial function: This sould not be used")
    # source and target must be of the same type
    T = typeof(source)
    for fld in fieldnames(source)
        if filedtype(T, fld) <: Array
            if isdefined(target,fld)
                getfield(target, fld)[:] = getfield(source, fld)[:]
            else
                setfield!(target, fld, copy(getfield(source, fld)))
            end
        else
            setfield!(target, fld, getfield(source, fld))
        end
    end
end


# Reads material parameters from a json file
function read_prms(filename::String)

    # read file
    file = open(filename, "r")
    data = JSON.parse(file)
    close(file)

    # parse materials
    mats_prms = Dict{String, Any}()
    for d in data
        name = d["name"]
        keys = d["prms"]
        vals = d["vals"]
        prms = Dict{Symbol, Float64}( Symbol(k) => v for (k,v) in zip(keys, vals) )
        mats_prms[name] = prms
    end

    return mats_prms
end


struct MaterialBind
    mat  :: Material
    expr :: Union{Symbol,Expr}
    function MaterialBind(target::Union{Symbol, Int, Array{Int,1}, UnitRange{Int}, Expr, String}, mat::Material)
        expr = :()
        ty = typeof(target)
        if ty==Symbol
            targets = (:all, :solids, :lines, :joints, :joints1D)
            !(target in targets) && error("MaterialBind: target should be one of $targets. :$target received")
            expr = target
        elseif ty==Expr
            expr = target
        #elseif ty==Int || ty==Array{Int,1} || ty==UnitRange{Int}
            #expr = :(id in $target)
        elseif ty<:Union{Int,String}
            #expr = :(tag==$target)
            expr = :(isequal(tag, $target))
        end

        return new(mat, expr)
    end
end
