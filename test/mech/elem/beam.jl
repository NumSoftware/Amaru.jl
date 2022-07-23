using Amaru
using Test

# 2D Truss

coord = [ 0 0; 1 0 ]
conn  = [ 1 2 ]

blt = BlockTruss(coord, conn, tag="beam")
msh = Mesh(blt, printlog=false)

# Finite element model

mats = [ "beam" => ElasticBeam(E=10, A=1, I=1) ]

bcs =
    [
     :(x==0 && y==0) => NodeBC(ux=0, uy=0, rz=0),
     :(x==1 && y==0) => NodeBC(fy=-10.),
    ]

#bc3 = BodyC( dom.elems, gz=-25 )
dom = Domain(msh, mats)

@test solve!(dom, bcs, printlog=false).success

