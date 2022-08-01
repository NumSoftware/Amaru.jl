using Amaru
using Test

printstyled("\nMesh revolve\n", color=:blue, bold=true)

shapes = (LIN2, LIN3)
data = ((48,49),(48,145)) 
for i in 1:length(shapes)
    shape = shapes[i]
    println("\nrevolving $(shape.name)")
    bl = Block( [0 0; 1 1], n=4, cellshape=shape, tag="solids")
    mesh = Mesh(bl, printlog=false)
    mesh = revolve(mesh, n=12)
    TR = @test (length(mesh.elems), length(mesh.nodes)) == data[i]
    println(TR)
end

shapes = (TRI3, TRI6, QUAD4, QUAD8)
data = ((768,437), (384,1113), (192,245), (192,921))
for i in 1:length(shapes)
    shape = shapes[i]
    println("\nrevolving $(shape.name)")
    bl = Block( [0 0; 1 1], nx=4, ny=4, cellshape=shape, tag="solids")
    mesh = Mesh(bl, printlog=false)
    mesh = revolve(mesh, n=12)
    TR = @test (length(mesh.elems), length(mesh.nodes)) == data[i]
    println(TR)
end
