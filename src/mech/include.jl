# This file is part of Amaru package. See copyright license in https://github.com/NumSoftware/Amaru

export elem_config_dofs, elem_init, elem_stiffness, elem_mass, elem_update!, elem_vals
export set_state
export solve!

include("elem/mech.jl")

# Elements
include("elem/mech-solid.jl")
#include("elem/in-mech-solid.jl")
include("elem/mech-rod.jl")
include("elem/mech-embrod.jl")
include("elem/mech-beam.jl")
include("elem/mech-joint.jl")
include("elem/mech-rsjoint.jl")
include("elem/mech-tipjoint.jl")
include("elem/mech-plate-mzc-quad4.jl")
include("elem/mech-plate-rm-quad4.jl")
include("elem/mech-plate-rm-quad8.jl")
include("elem/mech-shell-quad4.jl")
include("elem/mech-shell-degenerated.jl")

# Models for bulk elements
include("mat/elastic-solid.jl")
#include("mat/in-elastic-solid.jl")
include("mat/dp-solid.jl")
include("mat/vm-solid.jl")
include("mat/mazars-solid.jl")
include("mat/sc-solid.jl")
include("mat/compressiveconcrete-solid.jl")
include("mat/damageconcrete-solid.jl")
#include("mat/ortho-solid.jl")
#include("mat/ortho-concrete-solid.jl")

# Spring, dumper, lumped mass
include("elem/mech-lumpedmass.jl")
include("mat/lumpedmass.jl")
include("elem/mech-spring.jl")
include("mat/elastic-spring.jl")

# Models for truss elements
include("mat/elastic-rod.jl")
include("mat/pp-rod.jl")

# Models for beams
include("mat/elastic-beam.jl")

# Models for plates
include("mat/elastic-plateMZC.jl")
include("mat/elastic-plateRM.jl")
include("mat/elastic-plateRM8node.jl")

# Models for shels
include("mat/elastic-shell-quad4.jl")
include("mat/elastic-shell-degenerated.jl")

# Models for joint elements
include("mat/elastic-joint.jl")
include("mat/mc-joint.jl")
include("mat/p-joint.jl")
include("mat/mmc-joint.jl")
include("mat/tc-joint.jl")

# Models for 1D joint elements
include("mat/elastic-rsjoint.jl")
include("mat/ceb-rsjoint.jl")
include("mat/cyclic-rsjoint.jl")
include("mat/tipjoint.jl")
include("mat/elastic-tipjoint.jl")

# Stress-update integrator
include("mat/integrator.jl")

# Solvers
include("solver.jl")
include("dyn-solver.jl")
include("modal-solver.jl")

function reset_displacements(dom::Domain)
    for n in dom.nodes
        for dof in n.dofs
            for key in (:ux, :uy: :uz)
                haskey!(dof.vals, key) && (dof.vals[key]=0.0)
            end
        end
    end
end