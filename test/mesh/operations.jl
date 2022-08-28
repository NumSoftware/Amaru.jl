
using Amaru
using Test

printstyled("\nMesh move\n", color=:blue, bold=true)
bl = Block( [0 0; 1 1], nx=10, ny=10, cellshape=TRI3)
mesh = Mesh(bl)
move!(bl, dx=1)
mesh = Mesh(mesh, bl)
save(mesh, "out.vtk")
TR = @test length(mesh.nodes) == 431
println(TR)

printstyled("\nMesh extrude Mesh\n", color=:blue, bold=true)
bl = Block( [0 0; 1 1], nx=3, ny=3, cellshape=QUAD4)
mesh = Mesh(bl)
mesh = extrude(mesh, length=4, n=10)
save(mesh, "out.vtk")
TR = @test length(mesh.elems) == 90
println(TR)

printstyled("\nMesh extrude Block\n", color=:blue, bold=true)
bl = Block( [0 0; 1 1], nx=3, ny=3, cellshape=TRI3)
ble = extrude(bl, length=4, n=10)
mesh = Mesh(ble)
save(mesh, "out.vtk")
TR = @test length(mesh.elems) == 90
println(TR)

printstyled("\nMesh rotate\n", color=:blue, bold=true)
bl = Block( [0 0; 1 1], nx=4, ny=4, cellshape=QUAD4)
rotate!(bl, base = [0.5, 0.5], axis=[1,1], angle=45)
mesh = Mesh(bl)
save(mesh, "out.vtk")
TR = @test length(mesh.elems) == 16
println(TR)

printstyled("\nMesh polar\n", color=:blue, bold=true)
bl = Block( [0 0; 1 1], nx=4, ny=4, cellshape=QUAD4)
bls = polar(bl, base = [0, 0],  n=4)
mesh = Mesh(bls)
save(mesh, "out.vtk")
TR = @test length(mesh.elems) == 64
println(TR)

rm("out.vtk")
