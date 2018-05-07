using Amaru

# Mesh generation

blocks = [
    Block2D( [0 0; 1 2], nx=1, ny=4),
]

msh = Mesh(blocks, verbose=true)

# Finite element analysis

# Analysis data
k    = 1.0E-5  # permeability
E    = 5000.0  # Young modulus
nu   = 0.25    # Poisson
gw   = 10.0    # water specific weight

materials = [
    MaterialBind(:solids, LinSeep(k=k, gw=gw) ),
]

logger = [
            GroupLogger(:nodes, :(x==0) ),
]

dom = Domain(msh, materials, logger)

#for node in dom.nodes
    #@show node.id
    #@show node.dofs
#end

uw_f(t) = t/10.0
fw_f(t) = t/10.0

# Stage 1: loading

bcs = [
    #BC(:node, :(y==0), :(uw=$uw_f(t)) ),
    BC(:node, :(y==0), :(fw=$fw_f(t)) ),
    BC(:node, :(y==2), :(uw=0.) ),
]

hm_solve!(dom, bcs, end_time=500.0, tol=0.1, saveincs=true, verbose=true)

# Stage 2: draining
bcs = [
    #BC(:node, :(y==0), :(uw=$uw_f(t)) ),
    BC(:node, :(y==2), :(uw=0.) ),
]

hm_solve!(dom, bcs, end_time=1000.0, tol=0.1, nouts=1, saveincs=true, verbose=true)

#@show dom.nodes[44].dofs
#@show dom.nodes[1].dofs


# Output

using PyPlot
save(logger[1], "book.dat")

book = logger[1].book
@show length(book.tables)
@show book

for (i,table) in enumerate(book.tables)
    plot(table[:uw], table[:y], "-o")
end

@show book.tables[end][:uw]
#show()
