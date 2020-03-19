
using Amaru
using Test

printstyled("\nMesh move\n", color=:cyan)
bl = Block( [0 0; 1 1], nx=10, ny=10, cellshape=TRI3)
mesh = Mesh(bl, verbose=true)
move!(bl, dx=1)
mesh = Mesh(mesh, bl, verbose=true)
#save(mesh, "out.vtk")
@test length(mesh.points) == 231

printstyled("\nMesh extrude\n", color=:cyan)
bl = Block( [0 0; 1 1], nx=3, ny=3, cellshape=QUAD4)
mesh = Mesh(bl, verbose=true)
mesh = extrude(mesh, len=4, n=10)
#save(mesh, "out.vtk")
@test length(mesh.cells) == 90

bl = Block( [0 0; 1 1], nx=3, ny=3, cellshape=QUAD4)
ble = extrude(bl, len=4, n=10)
mesh = Mesh(ble, verbose=true)
#save(mesh, "out.vtk")
@test length(mesh.cells) == 90

printstyled("\nMesh rotate\n", color=:cyan)
bl = Block( [0 0; 1 1], nx=4, ny=4, cellshape=QUAD4)
rotate!(bl, base = [0.5, 0.5], axis=[1,1], angle=45)
mesh = Mesh(bl, verbose=true)
#save(mesh, "out.vtk")
@test length(mesh.cells) == 16

bl = Block( [0 0; 1 1], nx=4, ny=4, cellshape=QUAD4)
bls = polar(bl, base = [0, 0],  n=4)
mesh = Mesh(bls, verbose=true)
#save(mesh, "out.vtk")
@test length(mesh.cells) == 64

#rm("out.vtk")
