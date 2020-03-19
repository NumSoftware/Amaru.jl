using Amaru
using Test

# Mesh generation

blocks = [
    Block( [0 0; 1 2], nx=1, ny=4, tag="solids"),
]

msh = Mesh(blocks, silent=true)

# Finite element analysis

# Analysis data
k    = 100.62  # thermal conductivity w/m/k
rho  = 1.6   # water specific weight Ton/m3
cv   = 486e3  # specific heat (capacity) J/Ton/k

materials = [
             "solids" => LinThermo(k=k, rho=rho, cv=cv)
            ]
dom = Domain(msh, materials)

log1 = NodeGroupLogger()
loggers = [
    :(x==0) => log1
]
setloggers!(dom, loggers)


bcs = [
       :(y==0) => NodeBC(ut=10.0),
       :(y==2) => NodeBC(ut=20.0),
]

tm_solve!(dom, bcs, end_time=10000.0, tol=0.1, nincs=5, verbose=false)

# Output

if Amaru.config.makeplots
    using PyPlot
    save(log1, "book.dat")

    book = log1.book
    for (i,table) in enumerate(book.tables)
        plot(table[:ut], table[:y], "-o")
    end
end

