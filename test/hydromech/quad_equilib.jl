using Amaru

# Mesh generation

blocks = [
    Block2D( [0 0; 1 2], nx=10, ny=20),
]

msh = Mesh(blocks, verbose=true)

# Finite element analysis

# Analysis data
k    = 1.0E-5  # permeability
E    = 5000.0  # Young modulus
nu   = 0.25    # Poisson
gw   = 10.0    # water specific weight

materials = [
    MaterialBind(:solids, ElasticSolid(E=E, nu=nu) ),
]

logger = [
            GroupLogger(:nodes, :(x==0) ),
]

dom = Domain(msh, materials, logger)

#for node in dom.nodes
    #@show node.id
    #@show node.dofs
#end



t1 = 10.0

# Stage 1: loading

bcs = [
    BC(:node, :(x==0 && y==0), :(ux=0, uy=0) ),
    BC(:node, :(y==0), :(uy=0) ),
    BC(:face, :(y==2), :(ty=-10.0) ),
]

hm_solve!(dom, bcs, end_time=500.0, tol=0.1, saveincs=true, verbose=true)


# Output

using PyPlot
save(logger[1], "book.dat")

book = logger[1].book
@show length(book.tables)

for (i,table) in enumerate(book.tables)
    plot(table[:uy], table[:y], "-o")
end

@show book.tables[end][:uw]
#show()
