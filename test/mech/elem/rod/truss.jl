using Amaru
using Test


# 2D Truss

coord = [ 0 0; 9 0; 18 0; 0 9; 9 9; 18 9.]
conn  = [[1, 2], [1, 5], [2, 3], [2, 6], [2, 5], [2, 4], [3, 6], [3, 5], [4, 5], [5, 6]]

msh = Mesh(coord, conn)
tag!(msh.elems, "bars")

mats = [
        "bars" => ElasticRod(E=6.894757e7, A=0.043)
       ]

dom = Domain(msh, mats)

bcs = [
       :(x==0 && y==0) => NodeBC(ux=0, uy=0),
       :(x==0 && y==9) => NodeBC(ux=0, uy=0),
       :(x==9 && y==0) => NodeBC(fy=-450.),
       :(x==18&& y==0) => NodeBC(fy=-450.),
       :(x>=0)         => ElemBC(ty=-10.0)
      ]

@test solve!(dom, bcs, verbosity=1).success


# 3D Truss plus self weight

coord = [ 0.0 0.0 0.0; 0.0 1.0 0.0; 0.0 1.0 1.0]   # matriz de coordenadas
conn  = [[1, 3], [1, 2], [2, 3]]  # matriz de conectividades

msh = Mesh(coord, conn)
tag!(msh.elems, "bars")

dom = Domain(msh, mats)

bcs = [
       :(x==0 && y==0 && z==0) => NodeBC(ux=0),
       :(x==0 && y==1 && z==0) => NodeBC(ux=0, uy=0, uz=0),
       :(x==0 && y==1 && z==1) => NodeBC(ux=0, uy=0),
       :(x==0 && y==0) => NodeBC(fz=-50.),
       :(x>=0) => ElemBC(tz=-10.0)
      ]

@test solve!(dom, bcs, verbosity=0).success

save(dom, "dom.vtk")

