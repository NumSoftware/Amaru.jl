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
include("elem/mech-frame.jl")
include("elem/mech-solid.jl")
include("elem/mech-fluid.jl")
include("elem/mech-truss.jl")
include("elem/mech-embtruss.jl")
include("elem/mech-beam.jl")
include("elem/mech-joint.jl")
include("elem/mech-bondslip.jl")
include("elem/mech-sliptip.jl")
include("elem/mech-shell.jl")


# Models for bulk elements
include("mat/linear-elastic.jl")
include("mat/linear-elastic-fluid.jl")
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

# Models for joint and coohesive elements
include("mat/linear-joint.jl")
include("mat/linear-contact-joint.jl")
include("mat/mc-joint.jl")
include("mat/tc-cohesive.jl")
include("mat/asinh-yield-crack.jl")

# Models for 1D joint elements
include("mat/linear-bondslip.jl")
include("mat/ceb-bondslip.jl")
include("mat/cyclic-bondslip.jl")
include("mat/linear-sliptip.jl")
include("mat/linear-tipcontact.jl")

# Stress-update integrator
include("mat/integrator.jl")

# Solvers
include("mech-solver.jl")
include("dyn-solver.jl")
include("mech-modal-solver.jl")

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