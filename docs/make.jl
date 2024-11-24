# Inside make.jl
push!(LOAD_PATH,"../src/")
using Amaru, Documenter

root = joinpath(dirname(pathof(Amaru)), "..", "docs")

makedocs(
    # remotes = nothing,
    root = root,
    modules  = [Amaru],
    sitename = "Amaru",
    pagesonly = true,    # only listed pages are included
    checkdocs = :none,   # :missing, :all
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        # assets = ["assets/luxor-docs.css"],
        collapselevel=1,
    ),
    pages = [
        "Introduction" => "index.md",
        # "Tutorial" => "tutorial/tutorial.md",
        # "Mesh" =>  "meshing/meshing.md",
        # "Model model" =>  "modelling/modelling.md",
        "Mechanical analysis" =>  "modelling/mech/mech.md",
        "Examples" =>  [
            "2D truss" => "examples/truss.md",
            "Solid beam" => "examples/solid-beam.md",
        ]
        # "Plotting" =>  "plotting/plotting.md",
    ],
    doctest = false,
    repo = "https://github.com/NumSoftware/Amaru.jl",
    # edit_branch = "main",

    # repo = "github.com/JuliaGraphics/LuxorManual.git"
    # repo = "github.com/NumSoftware/Amaru.git"
)
        
# ENV["GITHUB_REPOSITORY"] = "NumSoftware/Amaru.jl"
# ENV["GITHUB_EVENT_NAME"] = "push"
# ENV["GITHUB_REF"]        = "main"
# ENV["GITHUB_ACTOR"]      = "usr"
# ENV["GITHUB_TOKEN"]      = raw"${{ secrets.GITHUB_TOKEN }}"

deploydocs(;
    devbranch = "main",
    target    = "build",
    branch    = "gh-pages",
    repo      = "github.com/NumSoftware/Amaru.jl.git",
)