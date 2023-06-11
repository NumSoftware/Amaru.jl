using Amaru
using Test

# 2D Truss

coord = [ 0 0; 1 0; 2 0 ]
conn  = [ [1, 2], [2, 3] ]
# coord = [ 0 0; 0.5 0; 1 0 ]
# conn  = [ [1, 2], [2, 3] ]

# blt = BlockTruss(coord, conn, tag="beam")
msh = Mesh(coord, conn, tag="beam")

# Finite element model

mats = [ "beam" => ElasticBeam(E=10, A=1, I=1) ]

bcs =
    [
     :(x==0 && y==0) => NodeBC(ux=0, uy=0),
     :(x==2 && y==0) => NodeBC(ux=0, uy=0),
     2 => NodeBC(mz=1.),
     2 => NodeBC(fy=-1.),
     2 => BodyC(ty=-12.),
    ]

#bc3 = BodyC( model.elems, gz=-25 )
model = FEModel(msh, mats)

@test solve!(model).success