using Amaru
using Test

# Mesh generation
bl  = Block3D( [0 0 0; 1.0 6.0 1.0], nx=1, ny=5, nz=3)
bl1 = BlockInset( [0.2 0.2 0.2; 0.2 5.8 0.2], curvetype="polyline", embedded=true)
bl2 = copy(bl1)
move!(bl2, dx=0.6)
bls = [ bl, bl1, bl2 ]

msh = Mesh(bls, verbose=true)
#msh = Mesh(bl, verbose=true)

#save(msh, "msh.vtk")


# FEM analysis
mats = [
        MaterialBind(:solids, ElasticSolid(E=1.e4, nu=0.25) ),
        #MaterialBind(:lines , ElasticRod(E=1.e8, A=0.005) ),
        MaterialBind(:lines , PPRod(E=1.e8, A=0.005, sig_y=500e3) ),
       ]


dom = Domain(msh, mats)


bcs = [
       NodeBC(:(y==0 && z==0), :(ux=0, uy=0, uz=0)),
       NodeBC(:(y==6 && z==0), :(ux=0, uy=0, uz=0)),
       FaceBC(:(z==1), :(tz=-1000)),
      ]

@test solve!(dom, bcs, nincs=20, verbose=true)
save(dom, "dom1.vtk")
