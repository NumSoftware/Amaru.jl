using Amaru

# Mesh generation
# ===============
# Mesh generation
blocks = [
    Block3D( [0 0 0; 0.2 2.0 0.2], nx=2, ny=12, nz=2, cellshape=HEX8, tag="solids"),
]

msh = Mesh(blocks, report=1);
mplot(msh, "mesh.pdf", field="elem-id")


# Finite element modeling
# =======================

materials = [
             "solids" => ElasticSolid(E=36e6, nu=0.2, rho=24.0),
            ]

# Finite element domain
domain = Model(msh, materials)

loggers = [
           :(x==0.2 && y==2 && z==0.2) => NodeLogger("end-node.dat"),
          ]

setloggers!(domain, loggers)

# List of boundary conditions
bcs = [
       :(y==0 && z==0)   => NodeBC(ux=0, uy=0, uz=0),
       :(y==2 && z==0)   => NodeBC(uz=0),
       :(y==1 && z==0.2) => NodeBC(fz=-10),
]

# Perform the finite element analysis
dynsolve!(domain, bcs, time_span=10.0, nincs=4, filekey="dyn")

save(domain, "domain.vtk")
mplot(domain, "domain.pdf", field="uy")
