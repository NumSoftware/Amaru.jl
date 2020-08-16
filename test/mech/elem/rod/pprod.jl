using Amaru
using Test

coord = [ 0. 0.; 1. 0. ]
conn  = [ [1, 2] ]

msh = Mesh(coord, conn, silent=true)
tag!(msh.elems, "bars")

mats = [
        "bars" =>  PPRod(E=210e6, A=0.01, sig_y=500e3),
        #"bars" => PPRod(E=210e6, A=0.01, sig_y=500e3, H=1000),
       ]

dom = Domain(msh, mats)

bcs = [
       :(x==0 && y==0) => NodeBC(ux=0, uy=0),
       :(x==1 && y==0) => NodeBC(uy=0),
       :(x==1 && y==0) => NodeBC(ux=0.003),
      ]


@test solve!(dom, bcs, nincs=10, verbose=false).success

