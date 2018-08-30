using Amaru
using Test

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
    NodeGroupLogger(:(x==0)),
]

dom = Domain(msh, materials, logger)

fw_f(t) = t/10.0

bcs = [
    NodeBC(:(y==0), :(fw=$fw_f(t)) ),
    NodeBC(:(y==2), :(uw=0.) ),
]

hm_solve!(dom, bcs, end_time=500.0, tol=0.1, verbose=true)

# Output

if Amaru.Debug.makeplots
    using PyPlot
    save(logger[1], "book.dat")

    book = logger[1].book
    for (i,table) in enumerate(book.tables)
        plot(table[:uw], table[:y], "-o")
    end
end

