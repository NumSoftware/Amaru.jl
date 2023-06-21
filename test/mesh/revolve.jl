using Amaru
using Test

printstyled("\nMesh revolve\n", color=:blue, bold=true)

shapes = (LIN2, LIN3)
data = ((48,49),(48,145)) 
for i in eachindex(shapes)
    shape = shapes[i]
    println("\nrevolving $(shape.name)")
    bl = Block( [0 0; 1 1], n=4, cellshape=shape, tag="solids")
    mesh = Mesh(bl)
    mesh = revolve(mesh, base=[0,0,0], axis=[0,1,0], n=12)
    TR = @test (length(mesh.elems), length(mesh.nodes)) == data[i]
    println(TR)
end

shapes = (TRI3, TRI6, QUAD4, QUAD8)
data = ((768,437), (384,1113), (192,245), (192,921))
for i in eachindex(shapes)
    shape = shapes[i]
    println("\nrevolving $(shape.name)")
    bl = Block( [0 0; 1 1], nx=4, ny=4, cellshape=shape, tag="solids")
    mesh = Mesh(bl)
    mesh = revolve(mesh, base=[0,0,0], axis=[0,1,0], n=12)
    TR = @test (length(mesh.elems), length(mesh.nodes)) == data[i]
    println(TR)
end

data = ((4,9), (64, 193))

println("\nrevolving node")
mesh = revolve(Node([1, 0, 0]), minangle=0, maxangle=90, base=[0,0,0], axis=[0,1,0], n=4)
TR = @test (length(mesh.elems), length(mesh.nodes)) == data[1]
println(TR)

println("\nrevolving cord")
mesh = revolve(mesh, base=[0,0,0], axis=[0,0,1], n=16)
TR = @test (length(mesh.elems), length(mesh.nodes)) == data[2]
println(TR)
