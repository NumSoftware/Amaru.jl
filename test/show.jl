using Amaru

# Finite element entities
bl  = Block2D( [0 0; 1 1], nx=4, ny=4, shape=QUAD9)
msh = Mesh(bl, verbose=true)


mat = MaterialBind(:solids, ElasticSolid(E=100.0, nu=0.2) )

mons = [
        Logger(:node, :(x==1 && y==1), "node.dat"),
        GroupLogger(:nodes, :(y==1), "nodes.dat"),
       ]

dom = Domain(msh, mat, mons)

bcs = [
        BC(:node, :(y==0), ux=0, uy=0),
        BC(:face, :(y==1), ty=2),
      ]


solve!(dom, bcs, nincs=3)

table = loadtable("node.dat")
book  = loadbook("nodes.dat")

save(table, "node.dat")
save(book, "nodes.dat")

rm("node.dat")
rm("nodes.dat")


println(bl)
println(msh)


println(mat)
println(mons[1])
println(mons)
println(bcs[1])
println(bcs)
println(dom.nodes[1].dofs[1])
println(dom.nodes[1].dofs)
println(dom.nodes[1])
println(dom.nodes)
println(dom.elems[1])
println(dom.elems)
println(dom)

println(table)
println(book)

#dump(STDOUT, dom, 1, "")
#dump(STDOUT, dom.elems[1:2], 3, "  ")
