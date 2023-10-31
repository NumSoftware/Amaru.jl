# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

__precompile__()

"""
**Amaru.jl**

Amaru module implements functions and types to perform finite element analyses.

**Important data types**

Node, Element, Model, Dof, Ip, NodeBC, SurfaceBC

"""
module Amaru

    using StatsBase, Statistics, LinearAlgebra, StaticArrays, SparseArrays, Arpack, Gmsh
    using Printf, DelimitedFiles, DataStructures, Glob, DocStringExtensions, Dates
    using Cairo, LaTeXStrings, MathTeXEngine

    import DataStructures: OrderedDict, OrderedSet

    # Tools module
    include("tools/include.jl")

    # generic exports
    export max, min, sort, reset, getindex, sort, copy!, show

    abstract type AbstractPoint end
    abstract type AbstractCell end
    abstract type AbstractBlock<:AbstractCell end
    abstract type AbstractDomain end


    # FEM
    include("dof.jl")
    export Dof, add_dof

    include("node.jl")
    export Node, get_data, setvalue!

    # Mesh
    include("mesh/include.jl")
    export getcoords

    include("analysis.jl")

    include("model-env.jl")
    export ModelEnv

    include("properties.jl")

    include("material.jl")
    export Material, read_prms

    include("ip.jl")
    export Ip, ip_vals, maximum, minimum, sort

    include("element.jl")
    export Element
    export getnodes, changequadrature!, get_ips, elems_ip_vals, update_material!, setstate!

    include("tag.jl")
    export tag!

    # Plotting
    include("plot/include.jl")
    export Chart, LinePlot, MeshPlot
    export addplot!

    # Boundary conditions
    include("bc.jl")
    export NodeBC, SurfaceBC, EdgeBC, BodyC

    include("logger.jl")
    export NodeLogger, IpLogger, SurfaceLogger
    export NodeSumLogger, FaceLogger, EdgeLogger
    export NodeGroupLogger, IpGroupLogger
    export PointLogger, SegmentLogger

    include("monitor.jl")
    export NodeMonitor
    export NodeSumMonitor
    export IpMonitor
    export IpGroupMonitor

    include("stage.jl")
    export Stage

    include("fe-model.jl")
    export Model
    export FEModel
    export Domain
    export addlogger!, addmonitor!, addloggers!, addmonitors!
    export setloggers!, setmonitors!

    include("solver.jl")

    include("io.jl")

    # Mechanical module
    include("mech/include.jl")

    # Hydromech module
    include("hydromech/include.jl")

    # Thermomech module
    include("thermomech/include.jl")

    # Acoustic module
    include("acousticmech/include.jl")

    include("deprecated.jl")

    # show function for FE related types
    Base.show(io::IO, obj::Dof)            = _show(io, obj, 2, "")
    Base.show(io::IO, obj::Node)           = _show(io, obj, 2, "")
    Base.show(io::IO, obj::Ip)             = _show(io, obj, 2, "")
    Base.show(io::IO, obj::IpState)        = _show(io, obj, 2, "")
    Base.show(io::IO, obj::Element)        = _show(io, obj, 2, "")
    Base.show(io::IO, obj::Material)       = _show(io, obj, 2, "")
    Base.show(io::IO, obj::BC)             = _show(io, obj, 2, "")
    Base.show(io::IO, obj::Facet)          = _show(io, obj, 2, "")
    Base.show(io::IO, obj::AbstractLogger) = _show(io, obj, 2, "")
    Base.show(io::IO, obj::Model)         = _show(io, obj, 2, "")

    # Chaining
    Base.:(<<)(a, b::Type{<:Element}) = return (a=>b)
    Base.:(<<)(a::Pair, b::Type{<:Material}) = return (a.first => a.second => b)
    Base.:(<<)(a::Pair, b::NamedTuple) = return (a.first => a.second.first => a.second.second => b)


    # testing
    export makeplots
    global makeplots=false
    export @runfiles
    macro runfiles(files)
        return esc(quote
            for file in $files
                printstyled("\nRunning file ", file,"...\n", color=:yellow, bold=true)
                include(file)
                println()
            end
        end)
    end

end#module
