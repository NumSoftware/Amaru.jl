using Amaru
using Test

printstyled("\nMesh generation on solids\n", color=:blue, bold=true)

println("\nGenerating mesh using TRI3")
bl = Block( [0 0; 1 1], nx=10, ny=10, cellshape=TRI3)
mesh = Mesh(bl, verbosity=0)
TR = @test length(mesh.nodes) == 121
println(TR)

println("\nGenerating mesh using TRI6")
bl = Block( [0 0; 1 1], nx=10, ny=10, cellshape=TRI6)
mesh = Mesh(bl, verbosity=0)
TR = @test length(mesh.nodes) == 441
println(TR)

println("\nGenerating mesh using QUAD4")
bl = Block( [0 0; 1 1], nx=10, ny=10, cellshape=QUAD4)
mesh = Mesh(bl, verbosity=0)
TR = @test length(mesh.nodes) == 121
println(TR)

println("\nGenerating mesh using QUAD8")
bl = Block( [0 0; 1 1], nx=10, ny=10, cellshape=QUAD8)
mesh = Mesh(bl, verbosity=0)
TR = @test length(mesh.nodes) == 341
println(TR)

println("\nGenerating mesh using QUAD9")
bl = Block( [0 0; 1 1], nx=10, ny=10, cellshape=QUAD9)
mesh = Mesh(bl, verbosity=0)
TR = @test length(mesh.nodes) == 441
println(TR)

println("\nGenerating mesh using HEX8")
bl = Block( [0 0 0; 1 1 1], nx=10, ny=10, nz=10, cellshape=HEX8)
mesh = Mesh(bl, verbosity=0)
TR = @test length(mesh.nodes) == 1331
println(TR)

println("\nGenerating mesh using HEX20")
bl = Block( [0 0 0; 1 1 1], nx=10, ny=10, nz=10, cellshape=HEX20)
mesh = Mesh(bl, verbosity=0)
TR = @test length(mesh.nodes) == 4961
println(TR)

println("\nGenerating mesh using HEX27")
bl = Block( [0 0 0; 1 1 1], nx=10, ny=10, nz=10, cellshape=HEX27)
mesh = Mesh(bl, verbosity=0)
TR = @test length(mesh.nodes) == 9261
println(TR)

println("\nGenerating mesh using TET4")
bl = Block( [0 0 0; 1 1 1], nx=10, ny=10, nz=10, cellshape=TET4)
mesh = Mesh(bl, verbosity=0)
TR = @test length(mesh.nodes) == 1331
println(TR)

println("\nGenerating mesh using TET10")
bl = Block( [0 0 0; 1 1 1], nx=10, ny=10, nz=10, cellshape=TET10)
mesh = Mesh(bl, verbosity=0)
TR = @test length(mesh.nodes) == 9261
println(TR)

# println("\nGenerating mesh using HEX8 in BlockCylinder")
# bl = BlockCylinder( [0 0 0; 5 5 5], r=2.0, nr=6, n=4, cellshape=HEX8)
# mesh = Mesh(bl, verbosity=0)
# TR = @test length(mesh.nodes) == 445
# println(TR)

# println("\nGenerating mesh using HEX20 in BlockCylinder")
# bl = BlockCylinder( [0 0 0; 5 5 5], r=2.0, nr=6, n=4, cellshape=HEX20)
# mesh = Mesh(bl, verbosity=0)
# TR = @test length(mesh.nodes) == 1641
# println(TR)

printstyled("\nMesh generation for surfaces\n", color=:blue, bold=true)

coords = [ 0.0 0.0 0.0
           1.0 0.0 0.0
           1.0 1.0 0.0
           0.0 1.0 0.0
           0.5 0.0 0.25
           1.0 0.5 0.0
           0.5 1.0 0.5
           0.0 0.5 0.0 ]
bl = Block(coords, nx=8, ny=8, cellshape=QUAD4)
mesh = Mesh(bl, verbosity=0)
TR = @test length(mesh.nodes) == 81

printstyled("\nMesh generation on trusses\n", color=:blue, bold=true)

coord = [ 0 0; 9 0; 18 0; 0 9; 9 9; 18 9.]
conn  = [ [1, 2], [1, 5], [2, 3], [2, 6], [2, 5], [2, 4], [3, 6], [3, 5], [4, 5], [5, 6] ]
mesh = Mesh(coord, conn)
TR = @test length(mesh.nodes) == 6
println(TR)

coord = [ 0.0 0.0 0.0; 0.0 1.0 0.0; 0.0 1.0 1.0]  
conn  = [ [1, 3], [1, 2], [2, 3] ]
mesh = Mesh(coord, conn)
TR = @test length(mesh.nodes) == 3
println(TR)

printstyled("\nMesh generation with insets\n", color=:blue, bold=true)

println("\nGenerating mesh using HEX8")
bl = Block( [0 0 0; 1 1 1], nx=8, ny=8, nz=8, cellshape=HEX8)
bli = BlockInset( [0 0 0; 1 1 1] )
mesh = Mesh(bl, bli)
TR = @test length(mesh.elems[:lines]) == 8
println(TR)

println("\nGenerating mesh using HEX20")
bl = Block( [0 0 0; 1 1 1], nx=8, ny=8, nz=8, cellshape=HEX20)
bli = BlockInset( [0 0 0; 1 1 1] )
mesh = Mesh(bl, bli)
TR = @test length(mesh.elems[:lines]) == 8
println(TR)

printstyled("\nMesh generation of joint cells\n", color=:blue, bold=true)

println("\nGenerating mesh using TRI3 and insets")
bl  = Block( [0 0; 1 1], nx=4, ny=4, cellshape=TRI3)
bli = BlockInset( [ 0 0; 1 1] )
mesh = Mesh(bl, bli, verbosity=0)
mesh = generate_joints!(mesh)
TR = @test length(mesh.elems) == 80
println(TR)

println("\nGenerating mesh using QUAD8 and insets")
bl  = Block( [0 0; 1 1], nx=4, ny=4, cellshape=QUAD8)
bli = BlockInset( [ 0 0; 1 1] )
mesh = Mesh(bl, bli, verbosity=0)
mesh = generate_joints!(mesh)
TR = @test length(mesh.nodes)==137 && length(mesh.elems)==48
println(TR)

println("\nGenerating mesh using QUAD8 and 3 layers")
bl  = Block( [0 0; 1 1], nx=4, ny=4, cellshape=QUAD8)
bli = BlockInset( [ 0 0; 1 1] )
mesh = Mesh(bl, bli, verbosity=0)
mesh = generate_joints!(mesh, layers=3)
@test length(mesh.nodes) == 158
TR = @test length(mesh.elems) == 48
println(TR)

println("\nGenerating mesh using TET4 and insets")
bl  = Block( [0 0 0; 1.0 2.0 1.0], nx=2, ny=4, nz=2, cellshape=TET4)
bli = BlockInset( [0 0 0; 1.0 2.0 1.0] )
mesh = Mesh(bl, bli, verbosity=0)
mesh = generate_joints!(mesh)
TR = @test length(mesh.elems) == 264
println(TR)
