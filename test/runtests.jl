# Constants
using Amaru
using Test

Amaru.config.makeplots = false

Base.show(io::IO, t::Test.Pass) = printstyled("\t[ ok ]", color=:green)
Base.show(io::IO, t::Test.Fail) = printstyled("\t[ fail ]", color=:red)

printstyled("\x1b[1m", "\nRunning tests...\n", "\x1b[0m", color=:green)
@testset begin
    include("mesh/runtests.jl")
    include("plot/runtests.jl")
    include("tools/runtests.jl")
    include("mech/runtests.jl")
    include("hydromech/runtests.jl")
    include("thermomech/runtests.jl")
end

include("clean.jl")
