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
    "mesh/plotting.jl",

    # Tools
    "tools/show.jl",
    "tools/logger.jl",

    # Static analysis: solid models
    "mech/solid/elastic-elems.jl",
    "mech/solid/elastic-hex8.jl",
    "mech/solid/dp.jl",
    "mech/solid/vm.jl",
    "mech/solid/mazars.jl",
    #"mech/solid/sc.jl",
    #"mech/solid/sc-beam.jl",

    # Static analysis: rod models
    "mech/rod/truss.jl",
    "mech/rod/pprod.jl",

    # Static analysis: beam models
    #"mech/beam/beam.jl",

    # Static analysis: joint models
    "mech/joint/joint1d.jl",
    "mech/joint/joint2d.jl",
    "mech/joint/joint2d_tag.jl",
    "mech/joint/joint3d_tag.jl",

    # Static analysis: embeddeed and semi-embedded elements
    "mech/inset/embedded.jl",
    "mech/inset/ceb.jl",

    # Dynamic analysis:
    "mech/dynamic/dyn-spring.jl",
    "mech/dynamic/dyn-solid.jl",

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
