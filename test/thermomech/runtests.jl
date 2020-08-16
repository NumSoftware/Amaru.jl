files = [    
         # Thermomechanical analysis: 
         "thermo.jl",
         "thermomech.jl",
        ]

for file in files
    printstyled("\nRunning file ", file,"...\n", color=:yellow, bold=true)
    include(file)
    println()
end
