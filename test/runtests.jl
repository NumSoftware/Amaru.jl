# Constants
using Amaru
using Test

path  = dirname(@__FILE__)
tests = readdir(path)
Amaru.config.makeplots = false

Base.show(io::IO, t::Test.Pass) = printstyled("\t[ ok ]", color=:green)
Base.show(io::IO, t::Test.Fail) = printstyled("\t[ fail ]", color=:red)

printstyled("\x1b[1m", "\nRunning tests...\n", "\x1b[0m", color=:green)

FILES = [
    # Meshing
    "mesh/shape_deriv.jl",
    "mesh/generation.jl",
    "mesh/io.jl",
    "mesh/operations.jl",
    "mesh/extrapolation.jl",
    "mesh/smoothing.jl",
    "mesh/revolve.jl",

    # Plotting
    "plot/mplot.jl",

    # Tools
    "tools/show.jl",
    "tools/logger.jl",

    # Static analysis: bulk elements
    "mech/elem/elastic-elems.jl",
    "mech/elem/elastic-hex8.jl",
    "mech/elem/axisymmetric.jl",
    "mech/mat/dp.jl",
    "mech/mat/vm.jl",
    "mech/mat/mazars.jl",
    #"mech/mat/sc.jl",
    #"mech/mat/sc-beam.jl",

    # Static analysis: rod models
    "mech/elem/rod/truss.jl",
    "mech/elem/rod/pprod.jl",

    # Static analysis: beam models
    #"mech/beam/beam.jl",

    # Static analysis: joint models
    "mech/elem/joint/joint1d.jl",
    "mech/elem/joint/joint2d.jl",
    #"mech/joint/joint2d_tag.jl",
    #"mech/joint/joint3d_tag.jl",

    # Static analysis: embeddeed and semi-embedded elements
    "mech/elem/inset/embedded.jl",
    "mech/elem/inset/ceb.jl",

    # Dynamic analysis:
    "mech/elem/dynamic/dyn-spring.jl",
    "mech/elem/dynamic/dyn-solid.jl",

    # Hydromechanical analysis:
    "hydromech/seep.jl",
    "hydromech/cutoff.jl",
    "hydromech/terzaghi.jl",
    "hydromech/terzaghi-joint.jl",
    "hydromech/drain.jl",
    "hydromech/drain-solid.jl",
    "hydromech/hm-drain.jl",

    # Thermomechanical analysis:
    "thermomech/thermo.jl",
    "thermomech/thermomech.jl",
   ]


@testset begin
    for file in FILES
        printstyled( "\nRunning file ", file,"...\n", color=:yellow, bold=true)
        #printstyled("Running test file ", file,"...\n", bold=true, color=:white)

        include(file)
        println()
    end
end

include("clean.jl")
