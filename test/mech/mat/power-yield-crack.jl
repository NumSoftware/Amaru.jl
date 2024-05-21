using Amaru
using Test

# Mesh
bl1  = Block( [0 0; 0.1 0; 0.1 0.1; 0 0.1], nx=1, ny=1, cellshape=QUAD4)
bl2  = Block( [0.1 0; 0.2 0; 0.2 0.1; 0.1 0.1], nx=1, ny=1, cellshape=QUAD4)
msh = Mesh(bl1, bl2)
insert_cohesive_elements!(msh, tag="joints")

tag!(msh.elems[BULKCELL, :(x<=0.125)].nodes, "left")
tag!(msh.elems[BULKCELL, :(x>=0.075)].nodes, "right")

# Finite element analysis
E  = 30.e6
ft = 3e3
fc = -30e3

mats = [
    :bulks => MechSolid => LinearElastic => (E=E, nu=0.2),
    :joints => MechJoint => PowerYieldCrack => (E=E, nu=0.2, fc=fc, ft=ft, zeta=5.0, wc=1.5e-4, alpha=1.5, gamma=0.1),
]

ana = MechAnalysis(stressmodel=:planestress, thickness=1.0)
model = FEModel(msh, mats, ana)

# Boundary conditions
bcs = [
    "left" => NodeBC(ux=0, uy=0),
    "right" => NodeBC(ux=-1e-9, uy=8e-5),
]
addstage!(model, bcs, nincs=20, nouts=20)

solve!(model, autoinc=true, maxits=3, tol=0.001, rspan=0.01, dTmax=0.1, scheme=:Ralston)
