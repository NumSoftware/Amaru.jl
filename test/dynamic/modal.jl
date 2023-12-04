using Amaru

geo = GeoModel()
s = 0.25
p1 = addpoint!(geo, 0, 0, 0, size=s)
p2 = addpoint!(geo, 1, 0, 0, size=s)
p3 = addpoint!(geo, 1, 10, 0, size=s)
p4 = addpoint!(geo, 0, 10, 0, size=s)

addline!(geo, p1, p2)
addline!(geo, p2, p3)
addline!(geo, p3, p4)
addline!(geo, p4, p1)

mesh = Mesh(geo)

# mplot(mesh, "mesh.pdf")


mats = [ :bulks => MechSolid => LinearElastic => (E=2e6, nu=0.2, rho=15.0) ]

ana = ModalAnalysis()

model = Model(mesh, mats, ana)

bcs = [
    :(y==0.0) => NodeBC(ux=0, uy=0)
]

addstage!(model, bcs)

solve!(model, quiet=false, nmods=5)

