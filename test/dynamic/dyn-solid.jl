using Amaru

# Mesh generation

blocks = [
    Block( [0 0 0; 0.2 2.0 0.2], nx=1, ny=10, nz=1, cellshape=HEX8, tag="solids"),
]

mesh = Mesh(blocks)

# Model definition

materials = [
    "solids" => MechSolid => LinearElastic => (E=30e6, nu=0.2, rho=24.0)
]

ctx = MechContext()
model = FEModel(mesh, materials, ctx)
ana = DynAnalysis(model)

# Finite element modeling
bcs = [
    :(y==0 && z==0)   => NodeBC(ux=0, uy=0, uz=0),
    :(y==2 && z==0)   => NodeBC(uz=0),
    :(y==1 && z==0.2) => NodeBC(fz=-10),
]

addstage!(ana, bcs, tspan=0.1, nincs=1, nouts=1)

solve!(ana, alpha=4.2038, beta=174.2803e-6, quiet=true)