using Amaru

geo = GeoModel()

p1 = Point(0, 0, 0, size=0.2)
p2 = Point(1, 0, 0, size=0.2)
p3 = Point(1, 1, 0, size=0.2)
p4 = Point(0, 1, 0, size=0.2)

addline!(geo, p1, p2)
addline!(geo, p2, p3)
addline!(geo, p3, p4)
addline!(geo, p4, p1)

pull!(geo, geo.surfaces, axis=[0,0,1], length=1)

# plotting
plot = GeometryPlot(geo, 
    lineweight = 0.5,
    fontsize = 8,
    # font = "Times New Roman",
)

save(plot, "geo.pdf")