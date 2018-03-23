# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

import Base.maximum
import Base.minimum
import Base.sort

"""
`IpState`

Abstract type for objects to store the state at integration points.
"""
abstract type IpState end



"""
`Ip(R, w)`

Creates an `Ip` object that represents an Integration Point in finite element analyses.
`R` is a vector with the integration point local coordinates and `w` is the corresponding integration weight.
"""
mutable struct Ip
    R    ::Array{Float64,1}
    w    ::Float64
    X    ::Array{Float64,1}
    id   ::Int
    tag  ::TagType
    owner::Any    # Element
    data ::IpState  # Ip current state
    data0::IpState  # Ip state for the last converged increment

    function Ip(R::Array, w::Float64)
        this     = new(vec(R), w)
        this.X   = Array{Float64}(3)
        this.tag = ""
        this.owner = nothing
        return this
    end
end

# The functions below can be used in conjuntion with sort
get_x(ip::Ip) = ip.X[1]
get_y(ip::Ip) = ip.X[2]
get_z(ip::Ip) = ip.X[3]


"""
`ip_vals`

Returns a dictionary with keys and vals for the integration point `ip`.
"""
function ip_vals(ip::Ip)
    coords = Dict( :x => ip.X[1], :y => ip.X[2], :z => ip.X[3] )
    vals   = ip_state_vals(ip.owner.mat, ip.data)
    return merge(coords, vals)
end

# Index operator for a ip collection using expression
function getindex(ips::Array{Ip,1}, cond::Expr) 
    condm = fix_comparison_scalar(cond)
    funex = :( (x,y,z,id,tag) -> false )
    funex.args[2].args[2] = condm
    fun = nothing
    try
        fun   = eval(funex)
    catch
        error("Ip getindex: Invalid condition ", cond)
    end

    result = Array{Ip}(0)
    for ip in ips
        x, y, z = ip.X
        if Base.invokelatest(fun, x, y, z, ip.id, ip.tag)
            push!(result, ip)
        end
    end

    length(result) == 0 && warn("getindex: No ips found for expression: $cond\n")
    return result
end

getindex(ips::Array{Ip,1}, cond::String) = getindex(ips, parse(cond))

# Get the maximum value of a given coordinate for the whole collection of ips
function maximum(ips::Array{Ip,1}, dir::Symbol) 
    idx = findfisrt((:x, :y, :z), dir)
    maximum([ip.X[idx] for ip in ips])
end

function minimum(ips::Array{Ip,1}, dir::Symbol) 
    idx = findfisrt((:x, :y, :z), dir)
    minimum([ip.X[idx] for ip in ips])
end

# Sort a collection of ips in a given direction
#function sort(ips::Array{Ip,1}, dir::Symbol=:x; rev::Bool=false) 
    #idx  = findfirst((:x, :y, :z), dir)
    #idxs = sortperm([ip.X[idx] for ip in ips], rev=rev)
    #return ips[idxs]
#end

