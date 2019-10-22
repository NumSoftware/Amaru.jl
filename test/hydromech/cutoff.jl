using Amaru

# Geometria 
Lrocha    = 30.00;
Hrocha    = 15.00;
Lgap      = 0.50;
Hgap      = 8.00;

blocks = [
    Block( [0 0; (Lrocha/2-Lgap/2) Hrocha], nx=14, ny=15, cellshape=QUAD8, tag="solids"),
    Block( [(Lrocha/2-Lgap/2) 0; (Lrocha/2+Lgap/2) (Hrocha-Hgap)], nx=1, ny=7, cellshape=QUAD8, tag="solids"),
    Block( [(Lrocha/2+Lgap/2) 0; Lrocha Hrocha], nx=14, ny=15, cellshape=QUAD8, tag="solids"),
]

mesh  = Mesh(blocks, silent=true) 

# Analysis data
k    = 1.0e-8  # permeability
gw   = 10.0    # water specific weight

tag!(mesh.cells[:solids], "solids")

# FEM analysis
mats = [
    "solids" => LinSeep(k=k, gammaw=gw, S=0.0),
]

dom = Domain(mesh, mats, gammaw=10.0)

# Boundary conditions
bcs = [
       :(x<=$Lrocha/2 && y==$Hrocha) => FaceBC(uw=10*gw),
       :(x>=$Lrocha/2 && y==$Hrocha) => FaceBC(uw=0),
      ]

# Solving
hm_solve!(dom, bcs, end_time=1, nincs=2, tol=1e-2, nouts=1, verbose=false)
