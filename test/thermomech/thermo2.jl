using Amaru
#using Test

# Mesh generation
blocks = [Block([0 0; 1 2], nx = 2, ny = 4, tag = "solids"),]

msh = Mesh(blocks,  report=false)

# Finite element analysis

# Analysis data
k = 50.2  # thermal conductivity W/m/K
rho = 7.8   # densidade ton/m3
cv = 486e3  # specific heat (capacity) J/ton/K

materials = ["solids" => LinThermo(k = k, rho = rho, cv = cv)]
model = Model(msh, materials)

log1 = NodeGroupLogger()
loggers = [:(y == 1) => log1]
setloggers!(model, loggers)


bcs = [:(x == 0) => NodeBC(ut = 10.0),

        :(y == 2) => NodeBC(ut = 20.0)

        ]

tm_solve!(model,end_time = 300000.0, tol = 0.1, nincs = 20)
# Output
save(model, "dom8.vtk")
save(log1, "book2.dat")
