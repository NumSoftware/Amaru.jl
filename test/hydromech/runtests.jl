files = [
         # Hydromechanical analysis
         "seep.jl",
         "cutoff.jl",
         "terzaghi.jl",
         "terzaghi-joint.jl",
         "drain.jl",
         "drain-solid.jl",
         "hm-drain.jl",
        ]

for file in files
    printstyled("\nRunning file ", file,"...\n", color=:yellow, bold=true)
    include(file)
    println()
end
