using Amaru
#using Test

# Mesh generation
blocks = [Block([0 0; 1 2], nx = 2, ny = 4, tag = "solids"),]

msh = Mesh(blocks,  verbosity=0)

# Finite element analysis

# Analysis data
k = 50.2  # thermal conductivity W/m/K
rho = 7.8   # densidade ton/m3
cv = 486e3  # specific heat (capacity) J/ton/K

materials = ["solids" => LinThermo(k = k, rho = rho, cv = cv)]
dom = Domain(msh, materials)

log1 = NodeGroupLogger()
loggers = [:(y == 1) => log1]
setloggers!(dom, loggers)


bcs = [:(x == 0) => NodeBC(ut = 10.0),

        :(y == 2) => NodeBC(ut = 20.0)

        ]

tm_solve!(dom, bcs, end_time = 300000.0, tol = 0.1, nincs = 20, verbosity=0)
# Output
save(dom, "dom8.vtk")
save(log1, "book2.dat")
