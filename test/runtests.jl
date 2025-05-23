using Amaru
using Test

# Base.show(io::IO, t::Test.Pass) = printstyled("\t[ ok ]", color=:green)
# Base.show(io::IO, t::Test.Fail) = printstyled("\t[ fail ]", color=:red)

printstyled("\x1b[1m", "\nRunning tests...\n", "\x1b[0m", color=:green)
Amaru.makeplots = false

let
    @testset begin
        @runfiles [
            "mesh/runtests.jl",
            "model/runtests.jl",
            "plot/runtests.jl",
            "tools/runtests.jl",
            "mech/runtests.jl",
            "dynamic/runtests.jl",
            "thermomech/runtests.jl",
            "acoustic-mech/runtests.jl",
            "hydromech/runtests.jl",
        ]
    end
end

include("clean.jl")
