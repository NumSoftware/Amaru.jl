
using Amaru
using Test

printstyled("\nExtrude Block\n", color=:blue, bold=true)
bl = Block( [0 0; 1 1], nx=3, ny=3, cellshape=TRI3)
ble = extrude(bl, axis=[0,0,1], length=4, n=10, quiet=true)
mesh = Mesh(ble)
TR = @test length(mesh.elems) == 90
println(TR)

printstyled("\nExtrude Mesh\n", color=:blue, bold=true)
bl = Block( [0 0; 1 1], nx=3, ny=3, cellshape=QUAD4)
mesh = Mesh(bl)
mesh = extrude(mesh, axis=[0,0,1], length=4, n=10, quiet=true)
TR = @test length(mesh.elems) == 90
println(TR)

printstyled("\nExtrude normal to Mesh\n", color=:blue, bold=true)

bl  = Block( [1 0 0; 1 1 0], n=6, cellshape=LIN3, tag="shell")
mesh = Mesh(bl)
mesh = revolve(mesh, minangle=45, maxangle=135, n=6, axis=[0,-1,0], base=[0,0,0])
mesh = extrude(mesh, length=4.25, n=6, quiet=true)
TR = @test length(mesh.elems) == 216
println(TR)
