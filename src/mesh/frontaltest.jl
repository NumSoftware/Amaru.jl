
using FemMesh

# Test

#bl = TBlock([0 0; 0.3 0; 1 0; 1 1; 0 1; 0 0.3], ns=[20,16,20,20,20,20])
#bl = TBlock([0 0; 0.3 0; 1 0; 1 1; 0 1; 0 0.3], ns=[20,16,10,10,16,20])
#bl = TBlock([0 0; 1 0; 1 1; 0 1], ns=[30,30,30,30])
#bl = TBlock([0 0; 1 0; 1 1], ns=[20,40,38])


#bl = TBlock([0 0; 1 0; 0.7 0.7; 0.1 1; 0 1], [ [1 2], [2 4 3], [4 5], [5 1] ], h=0.1, points=[1 0 0.05], lines=Dict(1=>0.1, 2=>0.2))
#bl = TBlock([0 0; 1 0; 0.7 0.7; 0.1 1; 0 1], [ [1 2], [2 3 4], [4 5], [5 1] ], h=0.1, sizehint=[1 0 0.05])
bl = TBlock([0 0; 1 0; 0.7 0.7; 0.1 1; 0 1], [ [1 2], [2 3], [3 4], [4 5], [5 1] ], h=0.1, sizehint=[1 0 0.05; 0 1 0.02; 0.5 0.5 0.1])


m = Mesh(bl)
save(m, "mesh1.vtk")
mplot(m, "mesh1.pdf")


smooth!(m, maxit=4, fixed=false)
#laplacian_smooth!(m, maxit=4, fixed=true)
save(m, "mesh2.vtk")
mplot(m, "mesh2.pdf")



#points = [ Point(0.0,0.0), Point(1.0,-1.0), Point(1.0,1.0), Point(0.0,2.0) ]

#boundary = disc_poly(points, 0.03)
#for node in boundary
    #@show "---"
    #@show node
    #@show node.prev
    #@show node.next
#end
#exit()

#points = [ ]
#k=1
#for node in boundary
    #global k+=1
    #node.data.id = k
    #push!(points, node.data)
#end


#mesh = Mesh()

#mesh.cells, mesh.points  = frontal(boundary, 0.03)
#mesh.points = get_points(mesh.cells)
#@show length(mesh.cells)

#fixup!(mesh)

#@show mesh.points

#save(mesh, "mesh1.vtk")
#smooth!(mesh, maxit=6)
#save(mesh, "mesh2.vtk")

