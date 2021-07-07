
using Amaru
using Test

printstyled("\nMesh move\n", color=:blue, bold=true)
bl = Block( [0 0; 1 1], nx=10, ny=10, cellshape=TRI3)
mesh = Mesh(bl, verbosity=0)
move!(bl, dx=1)
mesh = Mesh(mesh, bl, verbosity=0)
save(mesh, "out.vtk", verbosity=0)
TR = @test length(mesh.nodes) == 231
println(TR)

printstyled("\nMesh extrude Mesh\n", color=:blue, bold=true)
bl = Block( [0 0; 1 1], nx=3, ny=3, cellshape=QUAD4)
mesh = Mesh(bl, verbosity=0)
mesh = extrude(mesh, length=4, n=10)
save(mesh, "out.vtk", verbosity=0)
TR = @test length(mesh.elems) == 90
println(TR)

printstyled("\nMesh extrude Block\n", color=:blue, bold=true)
bl = Block( [0 0; 1 1], nx=3, ny=3, cellshape=TRI3)
ble = extrude(bl, length=4, n=10)
mesh = Mesh(ble, verbosity=0)
save(mesh, "out.vtk", verbosity=0)
TR = @test length(mesh.elems) == 90
println(TR)

printstyled("\nMesh rotate\n", color=:blue, bold=true)
bl = Block( [0 0; 1 1], nx=4, ny=4, cellshape=QUAD4)
rotate!(bl, base = [0.5, 0.5], axis=[1,1], angle=45)
mesh = Mesh(bl, verbosity=0)
save(mesh, "out.vtk", verbosity=0)
TR = @test length(mesh.elems) == 16
println(TR)

printstyled("\nMesh polar\n", color=:blue, bold=true)
bl = Block( [0 0; 1 1], nx=4, ny=4, cellshape=QUAD4)
bls = polar(bl, base = [0, 0],  n=4)
mesh = Mesh(bls, verbosity=0)
save(mesh, "out.vtk", verbosity=0)
TR = @test length(mesh.elems) == 64
println(TR)

printstyled("\nMesh 2D split using joints\n", color=:blue, bold=true)
bl  = Block( [0 0; 1 1], nx=4, ny=4, cellshape=TRI3)
bli = BlockInset( [ 0 0; 1 1] )
mesh = Mesh(bl, bli, verbosity=0)
mesh = generate_joints!(mesh)
TR = @test length(mesh.elems) == 80
println(TR)

printstyled("\nMesh 3D split using joints\n", color=:blue, bold=true)
bl  = Block( [0 0 0; 1.0 2.0 1.0], nx=2, ny=4, nz=2, cellshape=TET4)
bli = BlockInset( [0 0 0; 1.0 2.0 1.0] )
mesh = Mesh(bl, bli, verbosity=0)
mesh = generate_joints!(mesh)
save(mesh, "out.vtk", verbosity=0)
mesh = Mesh("out.vtk", verbosity=0)
save(mesh, "out.vtk", verbosity=0)
TR = @test length(mesh.elems) == 264
println(TR)

rm("out.vtk")
