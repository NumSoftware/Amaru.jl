using Amaru, Test

# 2D unstructured mesh
geo = GeoModel()

s= 0.3
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

# add a line that divides a surface
p1 = addpoint!(geo, 0, 0.5, 0, size=s)
p2 = addpoint!(geo, 1, 0.5, 0, size=s)
addline!(geo, p1, p2)

# add another line that divides a surface
p1 = addpoint!(geo, 0, 0.75, 0, size=s)
p2 = addpoint!(geo, 1, 0.75, 0, size=s)
addline!(geo, p1, p2)

mesh1 = Mesh(geo)

# Block based mesh
bl = Block([0 0; 1 1], nx=4, ny=4)
mesh2 = Mesh(bl)

meshes = [ mesh1, mesh2 ]

# FE analysis
mats = [
    :bulks => MechSolid => LinearElastic => (E=1e2, nu=0.499)
]
ana = MechAnalysis()

results = []
for mesh in [mesh1, mesh2]
    model = FEModel(mesh, mats, ana)
    
    log = NodeLogger()
    addlogger!(model, :(x==0.5 && y==1) => log )
    
    bcs = [
        :(y==0) => NodeBC(uy=0)
        :(x==0 && y==0) => NodeBC(ux=0)
        :(y==1) => SurfaceBC(ty=-1)
    ]
    
    addstage!(model, bcs)
    solve!(model)

    push!(results, log.table.uy[end])

end    

println(@test results[1] ≈ results[2] atol=1e-7)