# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

__precompile__()

"""
**Amaru.jl**

Amaru module implements functions and types to perform finite element analyses.

**Important data types**

Node, Element, Domain, Dof, Ip, NodeBC, FaceBC

"""
module Amaru
using Printf, StatsBase, Statistics, LinearAlgebra, SparseArrays, DelimitedFiles, Arpack
using DataStructures
import DataStructures.OrderedDict, DataStructures.OrderedSet

# Debug
mutable struct ConfigFlags
    debug::Bool
    makeplots::Bool
end
const config = ConfigFlags(false, true)

# eye function
eye(n::Int64) = Array{Float64}(I,n,n)

# Tools module
include("tools/include.jl")

# generic exports
export max, min, sort, reset, getindex, sort, copy!, show


#include("abstract.jl")
abstract type AbstractPoint end
abstract type AbstractCell end
abstract type AbstractMesh end


# FEM
include("model-env.jl")
export ModelEnv

include("node.jl")
export Node, Dof, add_dof, get_data

# Mesh
include("mesh/include.jl")

include("ip.jl")
export Ip, ip_vals, maximum, minimum, sort

include("material.jl")
export Material, read_prms

include("element.jl")
export Element
export get_nodes, changequadrature!, get_ips, elems_ip_vals, changemat!, setstate!

include("tag.jl")
export tag!


include("plot/mplot.jl")
include("plot/cplot.jl")

# Boundary conditions
include("bc.jl")
export NodeBC, FaceBC, EdgeBC, ElemBC

include("logger.jl")
export NodeLogger, IpLogger
export NodeSumLogger, FaceLogger, EdgeLogger
export NodeGroupLogger, IpGroupLogger
export PointLogger, SegmentLogger
export update_logger!

include("domain.jl")
export Domain, SubDomain, reset!, setloggers!

include("io.jl")

# Mechanical module
include("mech/include.jl")

# Hydromechanical module
include("hydromech/include.jl")

# Thermomechanical module
include("thermomech/include.jl")

# show function for FE related types
Base.show(io::IO, obj::Dof) = custom_dump(io, obj, 2, "")
Base.show(io::IO, obj::Node) = custom_dump(io, obj, 2, "")
Base.show(io::IO, obj::Ip) = custom_dump(io, obj, 2, "")
Base.show(io::IO, obj::IpState) = custom_dump(io, obj, 2, "")
Base.show(io::IO, obj::Element) = custom_dump(io, obj, 2, "")
Base.show(io::IO, obj::Material) = custom_dump(io, obj, 2, "")
Base.show(io::IO, obj::BC) = custom_dump(io, obj, 2, "")
Base.show(io::IO, obj::Facet) = custom_dump(io, obj, 2, "")
Base.show(io::IO, obj::AbstractLogger) = custom_dump(io, obj, 2, "")
Base.show(io::IO, obj::Domain) = custom_dump(io, obj, 2, "")

end#module
