using Amaru
using Test

coord = [ 0. 0.; 1. 0. ]
conn  = [ [1, 2] ]

msh = Mesh(coord, conn)
tag!(msh.elems, "bars")

mats = [
        # "bars" => MechBar => PPBar => (E=210e6, fy=500e3, A=0.01),
        # "bars" => MechBar => VonMises => (E=210e6, fy=500e3, A=0.01),
        "bars" => MechBar => PPRod => (E=210e6, fy=500e3, A=0.01),
       ]

ana = MechAnalysis()
model = FEModel(msh, mats, ana)

bcs = [
    :(x==0 && y==0) => NodeBC(ux=0, uy=0),
    :(x==1 && y==0) => NodeBC(uy=0),
    :(x==1 && y==0) => NodeBC(ux=0.003),
    ]
addstage!(model, bcs, nincs=10)


@test solve!(model).success

