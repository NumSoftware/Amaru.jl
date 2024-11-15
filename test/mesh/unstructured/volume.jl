using Amaru, Test

# 2D unstructured mesh with a hole
geo = GeoModel()

s= 0.1


#####

p1 = addpoint!(geo, 0, 0, 0, size=s)
p2 = addpoint!(geo, 2, 0, 0, size=s)
p3 = addpoint!(geo, 2, 2, 0, size=s)
p4 = addpoint!(geo, 0, 2, 0, size=s)

addline!(geo, p1, p2)
addline!(geo, p2, p3)
addline!(geo, p3, p4)
addline!(geo, p4, p1)

pull!(geo, geo.surfaces, axis=[0,0,1], length=1.0)

p5 = addpoint!(geo, 1, 0, 1, size=s)
p6 = addpoint!(geo, 1, 1, 1, size=s)
p7 = addpoint!(geo, 0, 1, 1, size=s)

addline!(geo, p5, p6)
addline!(geo, p6, p7)

pull!(geo, geo.surfaces[end], axis=[0,0,1], length=-0.5)

mesh = Mesh(geo)
save(mesh, "mesh.vtu")