files = [
         # Meshing
         "shape_deriv.jl",
         "generation.jl",
         "io.jl",
         "operations.jl",
         "extrapolation.jl",
         "smoothing.jl",
         "revolve.jl",
        ]

for file in files
    printstyled("\nRunning file ", file,"...\n", color=:yellow, bold=true)
    include(file)
    println()
end
