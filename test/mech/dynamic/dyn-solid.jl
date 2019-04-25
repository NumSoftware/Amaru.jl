using Amaru

# Mesh generation

blocks = [
    Block3D( [0 0 0; 0.2 2.0 0.2], nx=2, ny=12, nz=2, cellshape=HEX8, tag="solids"),
]

mesh = Mesh(blocks, verbose=true)

# Domain definition

materials = [
    "solids" => ElasticSolid(E=30e6, nu=0.2, rho=24.0)
]

dom = Domain(mesh, materials)

# Finite element modeling

bcs = [
    :(y==0 && z==0)   => NodeBC(ux=0, uy=0, uz=0),
    :(y==2 && z==0)   => NodeBC(uz=0),
    :(y==1 && z==0.2) => NodeBC(fz=-10),
]


dynsolve!(dom, bcs, time_span=0.1, nincs=50, verbose=true, filekey="dyn", nouts=50, autoinc=true, alpha=4.2038, beta=174.2803e-6)
#dynsolve!(dom, bcs, time_span=0.1, nincs=50, verbose=true, filekey="dyn", nouts=100, autoinc=true)
