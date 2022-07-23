using Amaru

# Mesh generation
# ===============

blocks = [
    Block3D( [0 0 0; 1 1 1], nx=1, ny=1, nz=1, cellshape=HEX8, tag="solids"),
]

msh = Mesh(blocks, verbosity=1);

#mplot(msh, "mesh.pdf", field="cell-id")


# Finite element modeling
# =======================

materials = [
             "solids" => ElasticSolid(E=100.0, nu=0.2),
            ]

# Finite element domain
domain = Domain(msh, materials)

loggers = [
           :(z==1) => FaceLogger("top-face.dat"),
          ]

setloggers!(domain, loggers)

# List of boundary conditions
bcs = [
       :(z==0) => NodeBC(ux=0, uy=0, uz=0),
       :(z==1) => SurfaceBC(tz=:(-10*x)),   # triangular load
]

# Perform the finite element analysis
solve!(domain, bcs, nincs=4, nouts=1)

save(domain, "domain.vtu")
mplot(domain, "domain.pdf", field="uy")
