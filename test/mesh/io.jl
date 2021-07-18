using Amaru
using Test

printstyled("\nWriting and loading vtk format\n", color=:blue, bold=true)
bl = Block( [0 0 0; 1 1 1], nx=10, ny=10, nz=10, cellshape=HEX8)
m1= Mesh(bl, printlog=false)

save(m1, "out.vtk", printlog=false)
sleep(0.1)
m2 = Mesh("out.vtk", printlog=false)
t = length(m1.nodes)==length(m2.nodes) && 
    length(m1.elems)==length(m2.elems) && 
    keys(m1.node_data)==keys(m2.node_data) &&
    keys(m1.elem_data)==keys(m2.elem_data) 

TR = @test t
println(TR)

save(m1, "out.vtu", printlog=false)
sleep(0.1)
m2 = Mesh("out.vtu", printlog=false)
t = length(m1.nodes)==length(m2.nodes) && 
    length(m1.elems)==length(m2.elems) && 
    keys(m1.node_data)==keys(m2.node_data) &&
    keys(m1.elem_data)==keys(m2.elem_data) 

TR = @test t
println(TR)
