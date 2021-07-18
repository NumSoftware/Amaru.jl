using Amaru
using Test

# Mesh generation
bl  = Block( [0 0 0; 1.0 6.0 1.0], nx=4, ny=6, nz=2, tag="solids")
bl1 = BlockInset( [0.2 0.2 0.2; 0.2 5.8 0.2], curvetype="polyline", tag="bars", jointtag="joints")
bl2 = move!( copy(bl1), dx=0.6)
bl3 = move!( copy(bl1), dx=0.3)
bls = [bl, bl1, bl2, bl3 ]

mesh = Mesh(bls, printlog=false)

# FEM analysis

mats = [
    "solids" => ElasticSolid(E=1.e4, nu=0.),
    "joints" => ElasticJoint1D(ks=1.e5, kn=1.e5, A=0.005),
    "bars"   => ElasticRod(E=1.e8, A=0.005),
]

dom = Domain(mesh, mats)

bcs = [
       :(y==0 && z==0) => NodeBC(uy=0, uz=0),
       :(y==6 && z==0) => NodeBC(uz=0),
       :(z==1) => FaceBC(tz=-1000 ),
      ]

mon = NodeLogger()
loggers = [
           :(x==0.5 && y==3.0 && z==0.5) => mon
          ]

setloggers!(dom, loggers)

@test solve!(dom, bcs, nincs=20, printlog=false).success

save(dom, "dom1.vtk")
