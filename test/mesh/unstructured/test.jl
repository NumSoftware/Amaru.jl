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

# pull!(geo, geo.faces[1], axis=[0,1,0], length=0.5)


# # add a line that divides a surface
# p1 = addpoint!(geo, 0, 0, 0.5, size=s)
# p2 = addpoint!(geo, 1, 0, 0.5, size=s)
# addline!(geo, p1, p2)


p1 = addpoint!(geo, 0.25, 0, 0.25, size=s)
p2 = addpoint!(geo, 0.75, 0, 0.25, size=s)
p3 = addpoint!(geo, 0.75, 0, 0.75, size=s)
p4 = addpoint!(geo, 0.25, 0, 0.75, size=s)

addline!(geo, p1, p2)
addline!(geo, p2, p3)
addline!(geo, p3, p4)
addline!(geo, p4, p1)


pull!(geo, geo.faces[1], axis=[0,1,0], length=0.5)
pull!(geo, geo.faces[2], axis=[0,1,0], length=0.5)


plot = GeometryPlot(geo); save(plot, "geo.pdf")