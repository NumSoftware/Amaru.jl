using Amaru

# Mesh generation

blocks = [
    Block2D( [0 0; 1 2], nx=1, ny=2),
]

mesh = Mesh(blocks, verbose=true)

# Finite element analysis

# Analysis data
k    = 1.0E-5  # permeability
E    = 5000.0  # Young modulus
nu   = 0.25    # Poisson
gw   = 10.0    # water specific weight

materials = [
    MaterialBind(:solids, ElasticSolidLinSeep(E=E, nu=nu, k=k, gw=gw) ),
]

logger = [
            GroupLogger(:nodes, :(x==0) ),
]

dom = Domain(mesh, materials, logger)


t1 = 10.0

# Stage 1: loading

bcs = [
    BC(:node, :(y==0), :(ux=0, uy=0) ),
    BC(:face, :(y==2), :(ty=-10.0) ),
    BC(:node, :(y==2), :(uw=0.) ),
]

hm_solve!(dom, bcs, end_time=10.0, saveincs=true, verbose=true)

# Stage 2: draining
bcs = [
    BC(:node, :(y==0), :(ux=0, uy=0) ),
    BC(:face, :(y==2), :(ty=-10.0) ),
    BC(:node, :(y==2), :(uw=0.) ),
]

hm_solve!(dom, bcs, end_time=1000.0, nincs=1, tol=1, nouts=1, saveincs=true, verbose=true)

#@show dom.nodes[44].dofs
#@show dom.nodes[1].dofs


# Output

using PyPlot
save(logger[1], "data.dat")

book = logger[1].book
@show length(book.tables)

for (i,table) in enumerate(book.tables)
    plot(table[:uw], table[:y], "-o")
end

@show book.tables[end][:uw]
show()
