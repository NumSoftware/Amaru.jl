using Amaru
using Base.Test

# Mesh generation
bl  = Block3D( [0 0 0; 1.0 6.0 1.0], nx=3, ny=20, nz=3)
bl1 = BlockInset( [0.2 0.2 0.2; 0.2 5.8 0.2], curvetype="polyline")
bl2 = move!( copy(bl1), dx=0.6)
bl3 = move!( copy(bl1), dx=0.3)
bls = [bl, bl1, bl2, bl3 ]

mesh = Mesh(bls, verbose=true)

# FEM analysis

mats = [
    MaterialBind(:solids  , ElasticSolid(E=1.e4, nu=0.) ),
    MaterialBind(:joints1D, ElasticJoint1D(ks=1.e5, kn=1.e5, A=0.005) ),
    MaterialBind(:lines   , ElasticRod(E=1.e8, A=0.005) ),
]

dom = Domain(mesh, mats)

bcs = [
       NodeBC(:(y==0 && z==0), :(uy=0, uz=0)),
       NodeBC(:(y==6 && z==0), :(uz=0)),
       FaceBC(:(z==1), :(tz=-1000 )),
      ]

mon = NodeLogger(:(x==0.5 && y==1.0 && z==0.5) )
#set_logger(dom, mon)

@test solve!(dom, bcs, nincs=20, verbose=true)

@show dom.nodes[:(y==3.0 && z==0)][1].dofdict[:uz].vals
save(dom, "dom1.vtk")
