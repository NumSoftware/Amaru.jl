using Amaru

# Mesh generation

blocks = [
    Block3D( [0 0 0; 1 1 1], nx=4, ny=4, nz=8, shape=HEX8),
]

mesh = Mesh(blocks, verbose=true)


# Domain definition

materials = [
    MaterialBind(:solids, ElasticSolid(E=100.0, nu=0.2) ),
    #MaterialBind(:(x>2), ElasticSolid(E=100.0, nu=0.2) ),
    #MaterialBind(:(tag=="concreto"), ElasticSolid(E=100.0, nu=0.2) ),
]

loggers = [
    Logger(:node, 1)  # by id
    #Logger(:ip, :(x>0))  # by expression
    #GroupLogger(:nodes, "nodes")  # by tag
    GroupLogger(:ips, :(x>0))  # by expression
]

dom = Domain(mesh, materials, loggers)


# Finite element modeling

bcs = [
    #BC(:node, :(x==0 && z==0), :(ux=0) ),
    #BC(:node, :(y==0 && z==0), :(uy=0) ),
    #BC(:node, :(z==0), :(uz=0) ),
    BC(:node, :(z==0), :(ux=0, uy=0, uz=0) ),
    #BC(:face, :(x==0 || x==1), ux=0),
    BC(:face, :(z==1), :(tz=-10.) ),
    #BC(:element, :all, :(tz=-10.) ),
    #BC(:face, :(z==1), :(tz=-10.*x) ),
]

solve!(dom, bcs, nincs=10, verbose=true)

save(dom, "dom1.vtk")

#println("  uz = ", dom.nodes[:(z==1)][1].dofdict[:uz].U)
println("  uz = ", dom.nodes[:(z==1)][1].dofdict[:uz].vals)
