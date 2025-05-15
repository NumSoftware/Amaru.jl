using Amaru

geo = GeoModel()
block = Block([0 0; 2 0 ], nx=2)
addblock!(geo, block)
mesh = Mesh(geo, ndim=2)

mat = [
    :lines => MechFrame => LinearElastic => (E=10, A=1, I=1),
      ]

ctx = MechContext()
model = FEModel(mesh, mat, ctx)

ana = MechAnalysis(model)

bcs = [ 
    x==0 => NodeBC(ux=0, uy=0),
    x==2 => NodeBC(ux=0, uy=0),
    x==1 => NodeBC(mz=1),
    x==1 => NodeBC(fy=-1),
    x>=1 => BodyC(qy=-12),
]

addstage!(ana, bcs, nincs=1, nouts=1)
run!(ana, autoinc=false)

    