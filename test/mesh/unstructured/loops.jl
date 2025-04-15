using Amaru
using Test

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

# inner loops
p5 = addpoint!(geo, 0.5, 0, 0, size=s)
p6 = addpoint!(geo, 0.5, 1, 0, size=s)
addline!(geo, p5, p6)

# hole loop
pa = addpoint!(geo, 0.2, 0.2, 0, size=s)
pb = addpoint!(geo, 0.4, 0.2, 0, size=s)
pc = addpoint!(geo, 0.4, 0.4, 0, size=s)
pd = addpoint!(geo, 0.2, 0.4, 0, size=s)
addline!(geo, pa, pb)
addline!(geo, pb, pc)
addline!(geo, pc, pd)
addline!(geo, pd, pa)

# outer loop
pa = addpoint!(geo, 1, 0.9, 0, size=s)
pb = addpoint!(geo, 1.1,0.9, 0, size=s)
pc = addpoint!(geo, 1.1,1.1, 0, size=s)
pd = addpoint!(geo, 0.9,1.1, 0, size=s)
pe = addpoint!(geo, 0.9,1, 0, size=s)
addline!(geo, pa, pb)
addline!(geo, pb, pc)
addline!(geo, pc, pd)
addline!(geo, pd, pe)

# outer loop
pa = addpoint!(geo, 0.4, 1, 0, size=s)
pb = addpoint!(geo, 0.4, 1.1, 0, size=s)
pc = addpoint!(geo, 0.2, 1.1, 0, size=s)
pd = addpoint!(geo, 0.2, 1, 0, size=s)
addline!(geo, pa, pb)
addline!(geo, pb, pc)
addline!(geo, pc, pd)

println(@test length(geo.faces) == 5)

# plot = GeometryPlot(geo); save(plot, "geo.pdf")




