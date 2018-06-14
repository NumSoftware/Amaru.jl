using Amaru
using Base.Test


# 2D Truss

coord = [ 0 0; 1 0 ]
conn  = [ 1 2 ]

blt = BlockTruss(coord, conn)
msh = Mesh(blt, verbose=false)

# Finite element model

mat = MaterialBind(:all, ElasticBeam(E=10, A=1, I=1) )

bcs = 
    [ 
     NodeBC(:(x==0 && y==0), :(ux=0, uy=0, rz=0))
     NodeBC(:(x==1 && y==0), :(fy=-10.))
    ]

#bc3 = ElemBC( dom.elems, gz=-25 )
dom = Domain(msh, mat)

@test solve!(dom, bcs, verbose=true)

