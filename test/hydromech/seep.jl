using Amaru
using Test

# Mesh generation

blocks = [
    Block2D( [0 0; 1 2], nx=1, ny=4, tag="solids"),
]

msh = Mesh(blocks, verbose=true)

# Finite element analysis

# Analysis data
k    = 1.0E-5  # permeability
E    = 5000.0  # Young modulus
nu   = 0.25    # Poisson
gw   = 10.0    # water specific weight

materials = [
    "solids" => LinSeep(k=k, gw=gw)
]
dom = Domain(msh, materials)

log1 = NodeGroupLogger()
loggers = [
    :(x==0) => log1
]
setloggers!(dom, loggers)


bcs = [
       :(y==0) => NodeBC(fw=:(t/10.0)),
       :(y==2) => NodeBC(uw=0.),
]

hm_solve!(dom, bcs, end_time=500.0, tol=0.1, verbose=true)

# Output

if Amaru.Debug.makeplots
    using PyPlot
    save(log1, "book.dat")

    book = log1.book
    for (i,table) in enumerate(book.tables)
        plot(table[:uw], table[:y], "-o")
    end
end

