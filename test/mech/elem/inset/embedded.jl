using Amaru
using Test

# Mesh generation
bl  = Block( [0 0 0; 1.0 6.0 1.0], nx=1, ny=5, nz=3, tag="solid")
bl1 = BlockInset( [0.2 0.2 0.2; 0.2 5.8 0.2], curvetype="polyline", embedded=true, tag="embedded")
bl2 = copy(bl1)
move!(bl2, dx=0.6)
bls = [ bl, bl1, bl2 ]

msh = Mesh(bls, silent=true)
#msh = Mesh(bl, silent=true)

#save(msh, "msh.vtk")


# FEM analysis
mats = [
        "solid" => ElasticSolid(E=1.e4, nu=0.25),
        "embedded" => PPRod(E=1.e8, A=0.005, sig_y=500e3),
        #"embedded" => ElasticSolid(E=1.e8, A=0.005),
       ]


dom = Domain(msh, mats)


bcs = [
       :(y==0 && z==0) => NodeBC(ux=0, uy=0, uz=0),
       :(y==6 && z==0) => NodeBC(ux=0, uy=0, uz=0),
       :(z==1) => FaceBC(tz=-1000),
      ]

@test solve!(dom, bcs, nincs=20, verbose=false).success
save(dom, "dom1.vtk")
