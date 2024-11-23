using Amaru, Test

# 3D unstructured mesh
geo = GeoModel()

s= 0.5
p1 = addpoint!(geo, 0, 0, 0, size=s)
p2 = addpoint!(geo, 1, 0, 0, size=s)
p3 = addpoint!(geo, 1, 0, 1, size=s)
p4 = addpoint!(geo, 0, 0, 1, size=s)

addline!(geo, p1, p2)
addline!(geo, p2, p3)
addline!(geo, p3, p4)
addline!(geo, p4, p1)

# add a line that divides a surface
p1 = addpoint!(geo, 0, 0, 0.5, size=s)
p2 = addpoint!(geo, 1, 0, 0.5, size=s)
addline!(geo, p1, p2)

pull!(geo, geo.surfaces, axis=[0, 1, 0], length=1)

# add a point to be embedded in top surface
addpoint!(geo, [0.5, 0.5, 1], size=s) 
mesh1 = Mesh(geo)

# Block based mesh
bl = Block([0 0 0; 1 1 1], nx=31, ny=3, nz=3)
mesh2 = Mesh(bl)

# FE analysis
mats = [
    :bulks => MechSolid => LinearElastic => (E=1e2, nu=0.499)
]
ctx = MechContext()

# FE analysis
results = []
for mesh in [mesh1, mesh2]
    model = FEModel(mesh, mats, ctx)
    ana = MechAnalysis(model)
    
    log = NodeLogger()
    addlogger!(ana, :(x==1 && y==1 && z==1) => log )
    
    bcs = [
        :(z==0) => NodeBC(uz=0)
        :(x==0 && y==0 && z==0) => NodeBC(ux=0, uy=0)
        :(z==1) => SurfaceBC(tz=-1)
    ]
    
    addstage!(ana, bcs)
    solve!(ana)

    push!(results, log.table.uz[end])

end    

println(@test results[1] ≈ results[2] atol=1e-3)



