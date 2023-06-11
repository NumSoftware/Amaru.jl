using Amaru
using Test

# Mesh generation

blocks = [
    Block( [0 0; 1 2], nx=1, ny=4, tag="solids"),
]

msh = Mesh(blocks)

# Finite element analysis

# Analysis data
k    = 1.0E-5  # permeability
E    = 5000.0  # Young modulus
nu   = 0.25    # Poisson
gw   = 10.0    # water specific weight

materials = [
    "solids" => LinSeep(k=k, gammaw=gw)
]
model = FEModel(msh, materials, gammaw=10)

log1 = NodeGroupLogger()
loggers = [
    :(x==0) => log1
]
setloggers!(model, loggers)

bcs = [
       :(y==0) => NodeBC(fw=:(t/10.0)),
       :(y==2) => NodeBC(uw=0.),
]

addstage!(model, bcs, tspan=500)
hm_solve!(model, tol=0.1)

# Output
if @isdefined(makeplots) && makeplots
    using PyPlot
    save(log1, "book.dat")

    book = log1.book
    for (i,table) in enumerate(book.tables)
        plot(table[:uw], table[:y], "-o")
    end
end

