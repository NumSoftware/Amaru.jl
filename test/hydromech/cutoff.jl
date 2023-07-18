using Amaru

# Geometry
Lrock    = 30.00;
Hrock    = 15.00;
Lgap      = 0.50;
Hgap      = 8.00;

blocks = [
    Block( [0 0; (Lrock/2-Lgap/2) Hrock], nx=14, ny=15, cellshape=QUAD8, tag="solids"),
    Block( [(Lrock/2-Lgap/2) 0; (Lrock/2+Lgap/2) (Hrock-Hgap)], nx=1, ny=7, cellshape=QUAD8, tag="solids"),
    Block( [(Lrock/2+Lgap/2) 0; Lrock Hrock], nx=14, ny=15, cellshape=QUAD8, tag="solids"),
]

mesh  = Mesh(blocks)

# Analysis data
k    = 1.0e-8  # permeability
gw   = 10.0    # water specific weight

tag!(mesh.elems.solids, "solids")

# FEM analysis
mats = [
    "solids" << SeepSolid << LinSeep << (k=k, S=0.0),
]

ana = HydroAnalysis(gammaw=gw)

model = FEModel(mesh, mats, ana)

# Boundary conditions
bcs = [
       :(x<=$Lrock/2 && y==$Hrock) << SurfaceBC(uw=10*gw),
       :(x>=$Lrock/2 && y==$Hrock) << SurfaceBC(uw=0),
      ]
addstage!(model, bcs, tspan=1, nincs=2, nouts=1)

# Solving
solve!(model, tol=1e-2)
