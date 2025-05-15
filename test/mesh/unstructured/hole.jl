using Amaru, Test

# 2D unstructured mesh with a hole
geo = GeoModel()

s= 0.03
p1 = addpoint!(geo, 0, 0, 0, size=s)
p2 = addpoint!(geo, 1, 0, 0, size=s)
p3 = addpoint!(geo, 1, 1, 0, size=s)
p4 = addpoint!(geo, 0, 1, 0, size=s)

addline!(geo, p1, p2)
addline!(geo, p2, p3)
addline!(geo, p3, p4)
addline!(geo, p4, p1)

# add a point at the middle top
addpoint!(geo, 0.5, 1, 0, size=s)

# add a hole
p1 = addpoint!(geo, 0.25, 0.25, 0, size=s)
p2 = addpoint!(geo, 0.75, 0.25, 0, size=s)
p3 = addpoint!(geo, 0.75, 0.75, 0, size=s)
p4 = addpoint!(geo, 0.25, 0.75, 0, size=s)

addline!(geo, p1, p2)
addline!(geo, p2, p3)
addline!(geo, p3, p4)
addline!(geo, p4, p1)

face = pickface(geo, 0.5, 0.5)
delete!(geo, face)

plot = GeometryPlot(geo); save(plot, "geo.pdf")

# Block based mesh
mesh1 = Mesh(geo, quadratic=true)

bl = Block([0 0; 0.25 0.25], nx=3, ny=3, cellshape=QUAD8)
bl_corners = array(bl, nx=2, ny=2, dx=0.75, dy=0.75)
bl = Block([0.25 0; 0.75 0.25], nx=6, ny=3, cellshape=QUAD8)
bl_bu = array(bl, nx=1, ny=2, dy=0.75)
bl = Block([0 0.25; 0.25 0.75], nx=3, ny=6, cellshape=QUAD8)
bl_lr = array(bl, nx=2, ny=1, dx=0.75)

mesh2 = Mesh(bl_corners, bl_bu, bl_lr)

# FE analysis
meshes = [ mesh1, mesh2 ]

mats = [
    :bulks => MechSolid => LinearElastic => (E=1e2, nu=0.25)
]
ctx = MechContext()

results = []
for mesh in [mesh1, mesh2]
    model = FEModel(mesh, mats, ctx)
    ana = MechAnalysis(model)
    log = NodeLogger()
    addlogger!(ana, :(x==0.5 && y==1) => log )
    
    bcs = [
        :(y==0) => NodeBC(uy=0)
        :(x==0 && y==0) => NodeBC(ux=0)
        :(y==1) => SurfaceBC(ty=-1)
    ]
    
    addstage!(ana, bcs)
    solve!(ana)

    push!(results, log.table.uy[end])
end    

println(@test results[1] ≈ results[2] atol=1e-3)