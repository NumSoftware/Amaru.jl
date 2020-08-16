files = [
         # Static analysis: bulk elements
         "elem/elastic-elems.jl",
         "elem/elastic-hex8.jl",
         "elem/axisymmetric.jl",
         "mat/dp.jl",
         "mat/vm.jl",
         "mat/mazars.jl",
         #"mat/sc.jl",
         #"mat/sc-beam.jl",

         # Static analysis: rod models
         "elem/rod/truss.jl",
         "elem/rod/pprod.jl",

         # Static analysis: beam models
         #"beam/beam.jl",

         # Static analysis: joint models
         "elem/joint/joint1d.jl",
         "elem/joint/joint2d.jl",
         #"joint/joint2d_tag.jl",
         #"joint/joint3d_tag.jl",

         # Static analysis: embeddeed and semi-embedded elements
         "elem/inset/embedded.jl",
         "elem/inset/ceb.jl",

         # Dynamic analysis:
         "elem/dynamic/dyn-spring.jl",
         "elem/dynamic/dyn-solid.jl",
        ]

for file in files
    printstyled("\nRunning file ", file,"...\n", color=:yellow, bold=true)
    include(file)
    println()
end
