using Amaru

# Mesh generation

blocks = [
    Block3D( [0 0 0; 1 1 1], nx=4, ny=4, nz=8, shape=HEX8),
]

mesh = Mesh(blocks, verbose=true)

tag!(mesh.points[1], 1)       # node 1
tag!(mesh.faces[:(z==0)], 10) # top face
tag!(mesh.faces[:(z==1)], 20) # bottom face


# Domain definition

materials = [
    MaterialBind(0, ElasticSolid(E=100.0, nu=0.2) ),
    #MaterialBind(:(x>2), ElasticSolid(E=100.0, nu=0.2) ),
    #MaterialBind(:(tag=="concreto"), ElasticSolid(E=100.0, nu=0.2) ),
    #MaterialBind(100, ElasticSolid(E=100.0, nu=0.2) ),
]

loggers = [
    Logger(:node, 1)  # by tag
    #Logger(:ip, :(x>0))  # by expression
    #GroupLogger(:nodes, "nodes")  # by tag
    #GroupLogger(:ips, :(x>0))  # by expression
]

dom = Domain(mesh, materials, loggers)

#set_logger!(dom, loggers)


# Finite element modeling

bcs = [
    BC(:face, 10, :(ux=0, uy=0, uz=0) ),
    BC(:face, 20, :(tz=-10.) ),
    #BC(:element, :all, :(tz=-10.) ),
    #BC(:face, :(z==1), :(tz=-10.*x) ),
]


solve!(dom, bcs, nincs=10, nouts=10, verbose=true)

save(dom, "dom.vtk")
rm("dom.vtk")


#println("  uz = ", dom.nodes[:(z==1)][1].dofdict[:uz].U)
println("  uz = ", dom.nodes[:(z==1)][1].dofdict[:uz].vals)
