using Amaru

geo = GeoModel(size=1)
s = 0.2

p1 = addpoint!(geo, 0.0, 0.0, size=s)
p2 = addpoint!(geo, 1.0, 0.0, size=s)
p3 = addpoint!(geo, 1.0, 1.0, size=s)
p4 = addpoint!(geo, 0.0, 1.0, size=s)

addpath!(geo, :M, p1, :L, p2, :L, p3, :L, p4, :L, p1)

pb1 = Point(0.2, 0.2)
pb2 = Point(0.5, 0.5)
pb3 = Point(0.6, 0.5)
pb4 = Point(0.9, 0.5)
pb5 = Point(1.2, 0.5)
addsubpath!(geo, :M, pb1, :C, pb2, pb3, pb4, :L, pb5)

mesh = Mesh(geo)

save(mesh, "dowel.vtu")

# FEM
mats = [
        :bulks => MechSolid => LinearElastic => (E=24e2, nu=0.2),

        # :lines => MechTruss => LinearElastic => (E=200e6, A=0.00011),
        :lines => MechBeam => LinearElastic => (E=200e6, A=0.00011),

        :linejoints => MechLSJoint => LinearLSInterface => (kn=5000, ks=6000, p=0.25),
        # :linejoints => MechLSJoint => CebLSInterface => (taumax=12, taures=3, s1=0.001, s2=0.0011, s3=0.004, alpha=0.5, beta=0.5, ks=(12/0.001)*5, kn=5000, p=0.25)
       ]

ana = MechAnalysis()
model = FEModel(mesh, mats, ana)

tag!(model.elems.lines.nodes[x>=1], "outside")

bcs = [
    y==0 => NodeBC(ux=0, uy=0),
    # "outside" => NodeBC(uy=0, ux=0.01),
    and(x==1.2, y==0.5) => NodeBC(ux=0.01, uy=-0.01),
]

addstage!(model, bcs, nincs=10)
solve!(model, tol=0.01, maxits=3, autoincr=true)
