using Amaru

# geo = GeoModel()

# p1 = Point(0, 0, 0, size=0.5)
# p2 = Point(1, 0, 0, size=0.5)
# p3 = Point(1, 1, 0, size=0.5)
# p4 = Point(0, 1, 0, size=0.5)

# addline!(geo, p1, p2)
# addline!(geo, p2, p3)
# addline!(geo, p3, p4)
# addline!(geo, p4, p1)

# pull!(geo, geo.faces, axis=[0,0,1], length=1.0)

geo = GeoModel()
block = Block([0 0 0; 1 1 1], nx=2, ny=2, nz=4, shape=HEX8)
addblock!(geo, block)

mesh = Mesh(geo, quiet=true)

# ctx = MechContext(stressmodel=:planestress)
ctx = MechContext()
mats = [
    :bulks => MechSolid => LinearElastic => (E=1e4, nu=0.25)
]
    
model = FEModel(mesh, mats, ctx, quiet=false)
ana = MechAnalysis(model)

bcs = [
    z==0 => NodeBC(ux=0, uy=0, uz=0)
    z==1 => SurfaceBC(tz=-100)
    # :(x<0.5 && y==1) => SurfaceBC(qy=-100)
]
# bcs = [
#     :(y==0) => NodeBC(ux=0,uy=0)
#     :(x<0.5 && y==1) => SurfaceBC(qy=-100)
# ]
addstage!(ana, bcs)
solve!(ana, quiet=false)

cmap = Colormap(:coolwarm, rev=true)

# plotting
plot = MeshPlot(model, 
    azimut = 30,
    elevation = 30,
    field = "uz",
    colormap = cmap,
    # warp = 20,
    # label = L"u_z", 
    # fontsize = 8,
    # font = "Times New Roman",
    # colorbarscale = 0.8
)

save(plot, "mesh.pdf")