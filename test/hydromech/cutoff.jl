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
    "solids" => LinSeep(k=k, gammaw=gw, S=0.0),
]

model = Model(mesh, mats, gammaw=10.0)

# Boundary conditions
bcs = [
       :(x<=$Lrock/2 && y==$Hrock) => SurfaceBC(uw=10*gw),
       :(x>=$Lrock/2 && y==$Hrock) => SurfaceBC(uw=0),
      ]
addstage!(model, bcs, tspan=1, nincs=2, nouts=1)

# Solving
hm_solve!(model, tol=1e-2, report=true)
