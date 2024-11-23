using Amaru
using Test

# Mesh generation

blocks = [
    Block( [0 0; 1 2], nx=1, ny=4, tag="solids"),
]

msh = Mesh(blocks)
# Finite element analysis

# Analysis data
k    = 100.62  # thermal conductivity w/m/k
rho  = 1.6   # water specific weight Ton/m3
cv   = 486e3  # specific heat (capacity) J/Ton/k

materials = [
             "solids" => ThermoSolid => ConstConductivity => (k=k, rho=rho, cv=cv)
]

ctx   = ThermoContext(T0=0)
model = FEModel(msh, materials, ctx)
ana   = ThermoAnalysis(model)

log1 = NodeGroupLogger()
loggers = [
    :(x==0) => log1
]
addloggers!(ana, loggers)


bcs = [
       :(y==0) => NodeBC(ut=10.0),
       :(y==2) => NodeBC(ut=20.0),
]
addstage!(ana, bcs, tspan=10000, nincs=5, nouts=1)
solve!(ana, tol=0.1)

# Output
if @isdefined(makeplots) && makeplots
    using PyPlot

    book = log1.book
    for (i,table) in enumerate(log1.book.tables)
        plot(table[:ut], table[:y], "-o")
    end
    show()
end

