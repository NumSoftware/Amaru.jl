# Inside make.jl
push!(LOAD_PATH,"../src/")
using Amaru, Documenter

makedocs(
         sitename = "Amaru",
         modules  = [Amaru],
         pages=[
                "Home" => "index.md"
               ])
            
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
    push_preview = true,
)