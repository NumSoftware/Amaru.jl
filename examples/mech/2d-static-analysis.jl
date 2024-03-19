using Amaru

# Mesh generation
blocks = [
    Block( [0 0; 3 0.4], nx=20, ny=10, cellshape=QUAD8, tag="solids"),
]

msh = Mesh(blocks)

# Finite element modeling
materials = [
             "solids" => MechSolid => LinearElastic => (E=200e6, nu=0.2),
            ]

ana = MechAnalysis(stressmodel=:planestress)
model = FEModel(msh, materials, ana)

addlogger!(model, :(x==1.5 && y==0) => NodeLogger("one-node.dat"))
addlogger!(model, :(y<0.025) => IpGroupLogger("ip-list.dat"))
addmonitor!(model, :(x==3 && y==0.4) => NodeMonitor(:uy))

bcs = [
    :(x==0 && y==0) => NodeBC(ux=0, uy=0),
    :(x==3 && y==0) => NodeBC(uy=0),
    :(y==0.4)       => SurfaceBC(ty=:(-0.1*x)), # triangular load
]
addstage!(model, bcs, nincs=5, nouts=5)

solve!(model)

# Plotting
mchart = MeshChart(model, 
    field = :sxx,
    colormap = :coolwarm,
    label = L"\sigma_x",
    warp = 20
)

save(mchart, "mesh.pdf")