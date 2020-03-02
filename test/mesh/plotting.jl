# Loading the finite element package Amaru
using Amaru
using Test


printstyled("\nPlotting\n", color=:blue, bold=true)

# 2D

blocks = [
    Block( [0 0; 3 0.4], nx=30, ny=8, cellshape=QUAD8),
]

print("Blocks 2D")
mplot(blocks, "out.pdf")
println( @test isfile("out.pdf") )
rm("out.pdf")

msh = Mesh(blocks, verbose=false);

print("Mesh 2D ")
mplot(msh, "out.pdf", field="quality")
println( @test isfile("out.pdf") )
rm("out.pdf")

# 3D

blocks = [
    Block( [0 0 0; 0.2 2.0 0.2], nx=2, ny=12, nz=2, cellshape=HEX8),
]
print("Blocks 3D")
mplot(blocks, "out.pdf")
println( @test isfile("out.pdf") )
rm("out.pdf")

msh = Mesh(blocks, verbose=false);

print("Mesh 3D ")
mplot(msh, "out.pdf", field="quality")
println( @test isfile("out.pdf") )
rm("out.pdf")
