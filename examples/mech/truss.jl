using Amaru

# 2D Truss
coord = [ 0 0; 9 0; 18 0; 0 9; 9 9; 18 9.]
conn  = [[1, 2], [1, 5], [2, 3], [2, 6], [2, 5], [2, 4], [3, 6], [3, 5], [4, 5], [5, 6]]

msh = Mesh(coord, conn, tag="bars")

mats = [
        "bars" => ElasticRod(E=6.894757e7, A=0.043)
       ]

model = Model(msh, mats, report=0)

bcs = [
       :(x==0 && y==0) => NodeBC(ux=0, uy=0),
       :(x==0 && y==9) => NodeBC(ux=0, uy=0),
       :(x==9 && y==0) => NodeBC(fy=-450.),
       :(x==18&& y==0) => NodeBC(fy=-450.),
       :(x>=0)         => BodyC(ty=-10.0)
      ]

solve!(model, bcs)

mplot(model, "truss.pdf", field="fa", colorbarlabel="axial force", colormap="gist_rainbow",  warpscale=100)
