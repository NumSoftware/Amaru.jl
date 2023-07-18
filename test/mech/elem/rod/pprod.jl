using Amaru
using Test

coord = [ 0. 0.; 1. 0. ]
conn  = [ [1, 2] ]

msh = Mesh(coord, conn)
tag!(msh.elems, "bars")

mats = [
        "bars" << MechRod << PPRod << (E=210e6, fy=500e3, A=0.01),
        #"bars" << PPRod << (E=210e6, A=0.01, sig_y=500e3, H=1000),
       ]

ana = MechAnalysis()
model = FEModel(msh, mats, ana)

bcs = [
    :(x==0 && y==0) << NodeBC(ux=0, uy=0),
    :(x==1 && y==0) << NodeBC(uy=0),
    :(x==1 && y==0) << NodeBC(ux=0.003),
    ]
addstage!(model, bcs, nincs=10)


@test solve!(model).success

