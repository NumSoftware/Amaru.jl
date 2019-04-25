using Amaru

# Finite element entities
bl  = Block2D( [0 0; 1 1], nx=4, ny=4, cellshape=QUAD9, tag="solids")
msh = Mesh(bl, verbose=true)


mats = [ "solids" =>  ElasticSolid(E=100.0, nu=0.2) ]
dom = Domain(msh, mats)

loggers = [
        :(x==1 && y==1) => NodeLogger("node.dat"),
        :(y==1)         => NodeGroupLogger("nodes.dat"),
       ]
setloggers!(dom, loggers)


bcs = [
        :(y==0) => NodeBC(ux=0, uy=0),
        :(y==1) => FaceBC(ty=2),
      ]


solve!(dom, bcs, nincs=3)

table = loadtable("node.dat")
book  = loadbook("nodes.dat")

rm("node.dat")
rm("nodes.dat")

println(bl)
println(msh)


println(mats)
println(loggers[1])
println(loggers)
println(bcs[1])
println(bcs)
println(dom.nodes[1].dofs[1])
println(dom.nodes[1].dofs)
println(dom.nodes[1])
println(dom.nodes)
println(dom.elems[1])
println(dom.elems)
println(dom)

#println(table)
#println(book)

#dump(STDOUT, dom, 1, "")
#dump(STDOUT, dom.elems[1:2], 3, "  ")
