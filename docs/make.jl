# Inside make.jl
push!(LOAD_PATH,"../src/")
using Amaru, Documenter

makedocs(
         sitename = "Amaru",
         modules  = [Amaru],
         pages=[
                "Home" => "index.md"
               ])
               
deploydocs(;
    repo="github.com/NumSoftware/Amaru.jl.git",
)