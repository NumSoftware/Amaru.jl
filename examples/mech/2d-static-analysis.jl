using Amaru

# Mesh generation
# ===============

blocks = [
    Block( [0 0; 3 0.4], nx=30, ny=50, cellshape=QUAD8, tag="solids"),
    # Block( [0 0; 3 0.4], nx=2, ny=4, cellshape=QUAD4, tag="solids"),
]

msh = Mesh(blocks);
# mplot(msh, "mesh.pdf", field="elem-id")

# Finite element modeling

materials = [
             "solids" => ElasticSolid(E=100.0, nu=0.2),
            ]

# Finite element model
# =====================

model = Model(msh, materials)

addlogger!(model, :(x==1.5 && y==0) => NodeLogger("one-node.dat"))
addlogger!(model, :(y<0.025) => IpGroupLogger("ip-list.dat"))
addmonitor!(model, :(x==3 && y==0.4) => NodeMonitor(:uy))

bcs = [
    :(x==0 && y==0) => NodeBC(ux=0, uy=0),
    :(x==3 && y==0) => NodeBC(uy=0),
    :(y==0.4)       => SurfaceBC(ty=:(-0.1*x)), # triangular load
]
addstage!(model, bcs, nincs=51, nouts=5)

# Perform the finite element analysis
solve!(model)

# save(model, "model.vtu")
# mplot(model, "model.pdf", field="uy")
