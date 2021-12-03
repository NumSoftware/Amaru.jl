# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

# Pressure Elements
include("elem/press.jl")
include("elem/press-solid.jl")

# Pressure - Mechanic Elements
include("elem/pressmech.jl")
include("elem/pressmech-interf.jl")

include("pressmech-solver.jl")

export pm_solve!
