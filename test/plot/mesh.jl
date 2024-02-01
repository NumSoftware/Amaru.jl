using Amaru

geo = GeoModel()

p1 = Point(0, 0, 0, size=0.1)
p2 = Point(1, 0, 0, size=0.1)
p3 = Point(1, 1, 0, size=0.1)
p4 = Point(0, 1, 0, size=0.1)

addline!(geo, p1, p2)
addline!(geo, p2, p3)
addline!(geo, p3, p4)
addline!(geo, p4, p1)

mesh = Mesh(geo)

ana = MechAnalysis(stressmodel="plane-stress")
mats = [
    :bulks => MechSolid => LinearElastic => (E=1e4, nu=0.25)
]
model = FEModel(mesh, mats, ana)

bcs = [
    :(y==0) => NodeBC(ux=0,uy=0)
    :(x<0.5 && y==1) => SurfaceBC(qy=-100)
]
addstage!(model, bcs)
solve!(model)

# plotting
plot = MeshPlot(model, 
    field = "uy",
    colormap = :coolwarm,
    warp = 20,
    label = L"u_z", 
    fontsize = 8,
    font = "Times New Roman",
    colorbarscale = 0.8
)

save(plot, "mesh.pdf")