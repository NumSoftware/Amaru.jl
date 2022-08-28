using Amaru
using Test

coord = [ 0. 0.; 1. 0. ]
conn  = [ 1 2 ]

blt = BlockTruss(coord, conn)
msh = Mesh(blt, silent=false)

mats = [
        MaterialBind(:all, PPRod(E=210e6, A=0.01, sig_y=500e3)),
        #MaterialBind(:all, PPRod(E=210e6, A=0.01, sig_y=500e3, H=1000)),
       ]

model = Model(msh, mats)

bcs = [
       NodeBC(:(x==0 && y==0), :(ux=0, uy=0)),
       NodeBC(:(x==1 && y==0), :(uy=0)),
       NodeBC(:(x==1 && y==0), :(ux=0.003)),
      ]


@test solve!(model,nincs=10).success

