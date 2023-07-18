# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

include("elem/hydromech.jl")
export elem_config_dofs, elem_init, elem_stiffness, update_elem!, elem_vals
export set_state

# Hydromech Elements
include("elem/hydromech-solid.jl")
include("elem/hydromech-joint.jl")
include("elem/hydromech-joint2.jl")

# Seep Elements
include("elem/seep-solid.jl")
include("elem/drainpipe.jl")
include("elem/seep-joint1d.jl")
include("elem/hydro-joint.jl")

# Models for solid elements (1D, 2D and 3D)
include("mat/elastic-solid-lin-seep.jl")
include("mat/lin-seep.jl")
include("mat/lin-drainpipe.jl")

# Models for joint elements
include("mat/elastic-joint-seep.jl")
include("mat/elastic-joint-seep2.jl")
include("mat/mc-joint-seep.jl")
include("mat/mc-joint-seep2.jl")
include("mat/tc-joint-seep.jl")
include("mat/joint1d-lin-seep.jl")
include("mat/joint-lin-seep.jl")

include("hydromech-solver.jl")
# include("hydromech-solver-bkp.jl")

export solve!
