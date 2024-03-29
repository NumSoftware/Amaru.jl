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
gw   = 10.0    # water specific weight

mats = [
    "solids" => SeepSolid => ConstPermeability => (k=k,)
]

ana = HydroAnalysis(gammaw=10.0)
model = FEModel(msh, mats, ana)

log1 = NodeGroupLogger()
loggers = [
    :(x==0) => log1
]
setloggers!(model, loggers)

bcs = [
       :(y==0) => NodeBC(fw=:(t/10.0)),
       :(y==2) => NodeBC(uw=0.0),
]

addstage!(model, bcs, tspan=500, nincs=1)
solve!(model, tol=0.1)

# Output
if @isdefined(makeplots) && makeplots
    using PyPlot
    save(log1, "book.dat")

    book = log1.book
    for (i,table) in enumerate(book.tables)
        plot(table[:uw], table[:y], "-o")
    end
end

