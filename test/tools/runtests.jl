files = [
         # Tools
         "show.jl",
         "logger.jl",
        ]

for file in files
    printstyled("\nRunning file ", file,"...\n", color=:yellow, bold=true)
    include(file)
    println()
end
