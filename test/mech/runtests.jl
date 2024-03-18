using Amaru

@runfiles [
         # Static analysis: bulk elements
         "elem/elastic-elems.jl",
         "elem/elastic-hex8.jl",
         "elem/elastic-quad4.jl",
         "elem/axisymmetric.jl",
         
         # Static analysis: rod elems
         "elem/rod/truss.jl",
         "elem/rod/pprod.jl",
         
         # Static analysis: beam elems
         "elem/beam/beam.jl",
         
         # Static analysis: shell elems
         "elem/shell/shell.jl",
         
         # Static analysis: joint models
         "elem/joint/joint1d.jl",
         "elem/joint/joint2d.jl",
         #"joint/joint2d_tag.jl",
         #"joint/joint3d_tag.jl",
         
         # Static analysis: embeddeed and semi-embedded elements
         "elem/inset/embedded.jl",
         "elem/inset/ceb.jl",
         "elem/inset/tip.jl",
         
         # Material models
         "mat/dp.jl",
         "mat/vm-2d.jl",
         "mat/vm-3d.jl",
         "mat/vm-beam-shell.jl",
         "mat/vm-beam-2d.jl",
         "mat/vm-beam-3d.jl",
         "mat/mazars.jl",
         #"mat/sc.jl",

        ]

