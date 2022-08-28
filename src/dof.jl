# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

"""
`Dof()`

Creates an object that represents a Degree of Freedom in a finite element analysis.
`Node` objects include a field called `dofs` which is an array of `Dof` objects.
"""
mutable struct Dof
    name      ::Symbol # essential variable name
    natname   ::Symbol # natural value name
    eq_id     ::Int64  # number of equation in global system
    prescribed::Bool   # flag for prescribed dof
    vals      ::OrderedDict{Symbol,Float64}
    function Dof(name::Symbol, natname::Symbol)
        new(name, natname, 0, false, OrderedDict{Symbol,Float64}())
    end
end


function Base.copy(dof::Dof)
    newdof = Dof(dof.name, dof.natname)
    newdof.eq_id = dof.eq_id
    newdof.prescribed = dof.prescribed
    newdof.vals = copy(dof.vals)
    return newdof
end


function Base.getindex(dofs::Array{Dof,1}, s::Symbol)
    for dof in dofs
        dof.name == s && return dof
        dof.natname == s && return dof
    end
    error("getindex: Dof key $s not found.")
end


function Base.haskey(dofs::Array{Dof,1}, s::Symbol)
    for dof in dofs
        dof.name == s && return true
        dof.natname == s && return true
    end
    return false
end


const null_Dof = Dof(:null, :null)
@inline null(::Type{Dof}) = null_Dof