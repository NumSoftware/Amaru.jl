# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export elem_config_dofs, elem_init, elem_stiffness, elem_mass, update_elem!, elem_vals
export set_state

@enum(StressState,
    PLANESTRESS = 0,
    PLANESTRAIN = 1,
    AXISYMMETRIC = 2,
    THREEDIMENSIONAL = 3,
)


include("elem/mech.jl")
include("elem/distributed.jl")

# Elements
include("elem/mech-solid.jl")
include("elem/mech-truss.jl")
include("elem/mech-embtruss.jl")
include("elem/mech-beam.jl")
include("elem/mech-joint.jl")
include("elem/mech-lsinterface.jl")
include("elem/mech-tipinterface.jl")
include("elem/mech-shell.jl")


# Models for bulk elements
include("mat/linear-elastic.jl")
include("mat/drucker-prager.jl")
include("mat/von-mises.jl")
include("mat/mazars.jl")
include("mat/smeared-crack.jl")
include("mat/damageconcrete-solid.jl")
include("mat/willam-warnke.jl")
include("mat/CSCP.jl")

# Spring, dumper, lumped mass
include("elem/mech-lumpedmass.jl")
include("mat/lumpedmass.jl")
include("elem/mech-spring.jl")
include("mat/linear-spring.jl")

# Models for plates
# include("mat/elastic-plateMZC.jl")
# include("mat/elastic-plateRM.jl")
# include("mat/elastic-plateRM8node.jl")

# Models for shells
# include("mat/elastic-shell-quad4.jl")
# include("mat/elastic-shell-degenerated.jl")
# include("mat/elastic-cook-shell.jl")

# Models for joint elements
include("mat/linear-joint.jl")

# Models for cohesive elements
include("mat/linear-cohesive.jl")
include("mat/mc-cohesive.jl")
# include("mat/p-cohesive.jl")
# include("mat/mmc-cohesive.jl")
include("mat/tc-cohesive.jl")
include("mat/tc2-cohesive.jl") 

# Models for 1D joint elements
include("mat/linear-lsinterface.jl")
include("mat/ceb-lsinterface.jl")
include("mat/cyclic-lsinterface.jl")
include("mat/linear-tipinterface.jl")
include("mat/nl-tipinterface.jl")

# Stress-update integrator
include("mat/integrator.jl")

# Solvers
include("mech-solver.jl")
include("dyn-solver.jl")
include("modal-solver.jl")

export solve!

# function reset_displacements(model::Model)
#     for n in model.nodes
#         for dof in n.dofs
#             for key in (:ux, :uy: :uz)
#                 haskey!(dof.vals, key) && (dof.vals[key]=0.0)
#             end
#         end
#     end
# end