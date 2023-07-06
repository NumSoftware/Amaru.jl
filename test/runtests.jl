using Amaru
using Test

# Base.show(io::IO, t::Test.Pass) = printstyled("\t[ ok ]", color=:green)
# Base.show(io::IO, t::Test.Fail) = printstyled("\t[ fail ]", color=:red)

# Tests

printstyled("\x1b[1m", "\nRunning tests...\n", "\x1b[0m", color=:green)

# files = [
#     "mesh/runtests.jl",
#     "plot/runtests.jl",
#     "tools/runtests.jl",
#     "mech/runtests.jl",
#     "dynamic/runtests.jl",
#     "hydromech/runtests.jl",
#     "thermomech/runtests.jl"
# ]
makeplots = false


let
    @testset begin
        # for file in files
        #     include(file)
        # end
        @runfiles [
            "mesh/runtests.jl",
            "plot/runtests.jl",
            "tools/runtests.jl",
            "mech/runtests.jl",
            "dynamic/runtests.jl",
            "hydromech/runtests.jl",
            "thermomech/runtests.jl"
        ]
        # runfiles(files)
    end
end

include("clean.jl")
