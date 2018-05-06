using Amaru
using Base.Test

coord = [ 0. 0.; 1. 0. ]
conn  = [ 1 2 ]

blt = BlockTruss(coord, conn)
msh = Mesh(blt, verbose=false)

mats = [ 
        MaterialBind(:all, PPRod(E=210e6, A=0.01, sig_y=500e3)),
        #MaterialBind(:all, PPRod(E=210e6, A=0.01, sig_y=500e3, H=1000)),
       ]

dom = Domain(msh, mats)

bcs = [
       BC(:node, :(x==0 && y==0), :(ux=0, uy=0)),
       BC(:node, :(x==1 && y==0), :(uy=0)),
       BC(:node, :(x==1 && y==0), :(ux=0.003)),
       #BC(:node, :(x==1 && y==0), :(fx=5010.0)),
      ]


@test solve!(dom, bcs, nincs=10, verbose=true)

