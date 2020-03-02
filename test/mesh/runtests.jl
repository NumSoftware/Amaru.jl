using Amaru
using Test

Base.show(io::IO, t::Test.Pass) = printstyled("\t[ ok ]", color=:green)
Base.show(io::IO, t::Test.Fail) = printstyled("\t[ fail ]", color=:red)

FILES = [
    "shape_deriv.jl"
    "generation.jl"
    "operations.jl"
    "extrapolation.jl"
    "smoothing.jl"
    "plotting.jl"
]

@testset begin
    for f in FILES
        printstyled( "\nRunning file ", f,"...\n", color=:yellow, bold=true)
        include(f)
    end
end




