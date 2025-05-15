using Amaru

geo = GeoModel()
s = 0.9

# p1 = addpoint!(geo, 0.0, 0.0, 0.0, size=s)
# p2 = addpoint!(geo, 1.0, 0.0, 0.0, size=s)
# p3 = addpoint!(geo, 1.0, 1.0, 0.0, size=s)
# p4 = addpoint!(geo, 0.0, 1.0, 0.0, size=s)

# addpath!(geo, :M, p1, :L, p2, :L, p3, :L, p4, :L, p1)

block = Block([0.0  0.0; 1.0 1.0], nx=3, ny=3)
addblock!(geo, block)
pa = addpoint!(geo, 0.1, 0.5, 0.0, size=s)
pb = addpoint!(geo, 0.9, 0.5, 0.0, size=s)
addsubpath!(geo, :M, pa, :L, pb, embedded=true)

# face = pickface(geo, 0.1, 0.1, 0.0)
# pull!(geo, face, axis=[0,0,1], length=0.5)

mesh = Mesh(geo)
# plot = GeometryPlot(geo); save(plot, "geo.pdf")
plot = MeshPlot(mesh, nodelabels=true); save(plot, "mesh.pdf")


# insert_cohesive_elements!(mesh)
# add_boundary_interface_elements!(mesh, z==0, tag="interface", supporttag="suport")