using Amaru
using Test

printstyled("\nWriting and loading vtk format\n", color=:blue, bold=true)
bl = Block( [0 0 0; 1 1 1], nx=10, ny=10, nz=10, cellshape=HEX8)
m1= Mesh(bl, silent=true)

save(m1, "out.vtk", verbose=false)
m2 = Mesh("out.vtk", verbose=false)
t = length(m1.points)==length(m2.points) && 
    length(m1.cells)==length(m2.cells) && 
    keys(m1.point_data)==keys(m2.point_data) &&
    keys(m1.cell_data)==keys(m2.cell_data) 

TR = @test t
println(TR)

save(m1, "out.vtu", verbose=false)
m2 = Mesh("out.vtu", verbose=false)
t = length(m1.points)==length(m2.points) && 
    length(m1.cells)==length(m2.cells) && 
    keys(m1.point_data)==keys(m2.point_data) &&
    keys(m1.cell_data)==keys(m2.cell_data) 

TR = @test t
println(TR)
