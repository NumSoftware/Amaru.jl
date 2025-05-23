using Amaru
using Test

# Mesh generation
bl  = Block( [0 0 0; 1.0 6.0 1.0], nx=4, ny=6, nz=2, tag="solids")
bl1 = BlockInset( [0.2 0.2 0.2; 0.2 5.8 0.2], curvetype="polyline", tag="bars", jointtag="joints")
bl2 = move!( copy(bl1), dx=0.6)
bl3 = move!( copy(bl1), dx=0.3)
bls = [bl, bl1, bl2, bl3 ]

mesh = Mesh(bls)

# FEM analysis

mats = [
    "solids" => MechSolid => LinearElastic => (E=1.e4, nu=0.),
    "joints" => MechBondSlip => ElasticRSJoint => (ks=1.e5, kn=1.e5, p=0.25), 
    "bars"   => MechBar => LinearElastic => (E=1.e8, A=0.005),
]

ctx = MechContext()
model = FEModel(mesh, mats, ctx)
ana = MechAnalysis(model)

bcs = [
       :(y==0 && z==0) => NodeBC(uy=0, uz=0),
       :(y==6 && z==0) => NodeBC(uz=0),
       :(z==1) => SurfaceBC(tz=-1000 ),
      ]

addlogger!(ana, :(x==0.5 && y==3.0 && z==0.5) => NodeLogger() )
addstage!(ana, bcs, nincs=20)
@test solve!(ana).success