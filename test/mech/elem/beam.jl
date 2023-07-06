using Amaru
using Test

# 2D Truss

coord = [ 0 0; 1 0 ]
conn  = [ 1 2 ]

blt = BlockTruss(coord, conn, tag="beam")
msh = Mesh(blt)

# Finite element model

mats = [ "beam" << ElasticBeam(E=10, A=1, I=1) ]

bcs =
    [
     :(x==0 && y==0) << NodeBC(ux=0, uy=0, rz=0),
     :(x==1 && y==0) << NodeBC(fy=-10.),
    ]

#bc3 = BodyC( model.elems, gz=-25 )
ana = MechAnalysis()
model = FEModel(msh, mats, ana)

@test solve!(model).success

