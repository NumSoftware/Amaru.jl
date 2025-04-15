using Amaru, Test

# 3D unstructured mesh
geo = GeoModel()

s= 0.5
p0 = addpoint!(geo, 0, 0, 0, size=s)
p1 = addpoint!(geo, 1, 0, 0, size=s)
p2 = addpoint!(geo, 0, 1, 0, size=s)
p3 = addpoint!(geo, -1, 0, 0, size=s)

addarc!(geo, p1, p0, p2)
addarc!(geo, p2, p0, p3)
addline!(geo, p1, p3)

pull!(geo, geo.faces[1], axis=[0,0,1], length=0.5)

# plot = GeometryPlot(geo); save(plot, "geo.pdf")