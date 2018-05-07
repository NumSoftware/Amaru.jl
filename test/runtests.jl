# Constants
using Amaru
using Base.Test

path  = dirname(@__FILE__)
tests = readdir(path)
const NOPLOTS = 1

print_with_color(:green, "\x1b[1m", "\nRunning tests...\n", "\x1b[0m")

FILES = [
    "fem/show.jl",
    "fem/logger.jl",

    # Static analysis: solid elements
    "mech/elastic-solid/solid.jl",
    "mech/elastic-solid/elastic.jl",
    "mech/elastic-solid/one-elem.jl",

    # Static analysis: joint elements
    "mech/joint/joint1d.jl",
    "mech/joint/joint2d.jl",

    # Static analysis: embeddeed and semi-embedded elements
    "mech/semi-embedded/embed.jl",
    "mech/semi-embedded/ceb.jl",

    # Nonlinear static analysis: solid elements
    "mech/nl-solid/dp.jl",
    "mech/nl-solid/vm.jl",
    "mech/nl-solid/mazars.jl",
    #"mech/nl-solid/sc.jl",
    #"mech/nl-solid/sc-beam.jl",

    # Dynamic analysis: solid elements
    #"mech/daynamic/dyn-solid.jl",

    # Hydromechanical analysis: solid elements
    #"hydromech/terzaghi.jl",
   ]


@testset begin
    for file in FILES
        print_with_color(:white, "Running test file ", file,"...\n", bold=true)
        include(file)
        println()
    end
end
