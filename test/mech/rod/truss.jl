using Amaru
using Test


# 2D Truss

coord = [ 0 0; 9 0; 18 0; 0 9; 9 9; 18 9.]
conn  = [ 1 2; 1 5; 2 3; 2 6; 2 5; 2 4; 3 6; 3 5; 4 5; 5 6]

blt = BlockTruss(coord, conn)
msh = Mesh(blt, verbose=false)

mats = [
        MaterialBind(:all, ElasticRod(E=6.894757e7, A=0.043)),
        #MaterialBind(:all, ElasticRod(E=6.894757e7, A=1.0)),
       ]

dom = Domain(msh, mats)

bcs = [
       NodeBC(:(x==0 && y==0), :(ux=0, uy=0)),
       NodeBC(:(x==0 && y==9), :(ux=0, uy=0)),
       NodeBC(:(x==9 && y==0), :(fy=-450.)),
       NodeBC(:(x==18&& y==0), :(fy=-450.)),
       ElemBC(:(x>=0), :(ty=-10.0))
      ]

@test solve!(dom, bcs, verbose=true)


# 3D Truss plus self weight

coord = [ 0.0 0.0 0.0; 0.0 1.0 0.0; 0.0 1.0 1.0]   # matriz de coordenadas
conn  = [ 1 3; 1 2;2 3]  # matriz de conectividades

blt = BlockTruss(coord, conn)
msh = Mesh(blt, verbose=false)

dom = Domain(msh, mats)

bcs = [
       NodeBC(:(x==0 && y==0 && z==0), :(ux=0)),
       NodeBC(:(x==0 && y==1 && z==0), :(ux=0, uy=0, uz=0)),
       NodeBC(:(x==0 && y==1 && z==1), :(ux=0, uy=0)),
       NodeBC(:(x==0 && y==0), :(fz=-50.)),
       ElemBC(:(x>=0), :(tz=-10.0))
      ]

@test solve!(dom, bcs, verbose=true)

save(dom, "dom.vtk")

