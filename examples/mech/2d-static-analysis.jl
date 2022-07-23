using Amaru

# Mesh generation
# ===============

blocks = [
    Block( [0 0; 3 0.4], nx=30, ny=8, cellshape=QUAD8, tag="solids"),
]

msh = Mesh(blocks, verbosity=1);

mplot(msh, "mesh.pdf", field="elem-id")

# Finite element modeling

materials = [
             "solids" => ElasticSolid(E=100.0, nu=0.2),
            ]

# Finite element domain
# =====================

domain = Domain(msh, materials)

loggers = [
           :(x==1.5 && y==0) => NodeLogger("one-node.dat"),
           :(y<0.025) => IpGroupLogger("ip-list.dat"),
          ]

setloggers!(domain, loggers)

# List of boundary conditions
bcs = [
       :(x==0 && y==0) => NodeBC(ux=0, uy=0),
       :(x==3 && y==0) => NodeBC(uy=0),
       :(y==0.4)       => SurfaceBC(ty=:(-0.1*x)), # triangular load
]

# Perform the finite element analysis
solve!(domain, bcs, nincs=10)

save(domain, "domain.vtu")
mplot(domain, "domain.pdf", field="uy")
