# Constants
using Amaru
using Test

path  = dirname(@__FILE__)
tests = readdir(path)
Amaru.debug.makeplots = false


printstyled("\x1b[1m", "\nRunning tests...\n", "\x1b[0m", color=:green)

FILES = [
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
    "hydromech/terzaghi.jl",
    "hydromech/terzaghi-joint.jl",
   ]


@testset begin
    for file in FILES
        printstyled("Running test file ", file,"...\n", bold=true, color=:white)
        @show file
        include(file)
        println()
    end
end

include("clean.jl")
