using Amaru
#using Test

# Mesh generation
blocks = [Block([0 0; 1 2], nx = 4, ny = 8, tag = "solids"),]

msh = Mesh(blocks)

# Finite element analysis

# Analysis data
k = 50.2  # thermal conductivity W/m/K
rho = 7.8   # densidade ton/m3
cv = 486e3  # specific heat (capacity) J/ton/K

materials = ["solids" << ConstConductivity(k = k, rho = rho, cv = cv)]
model = FEModel(msh, materials)

log1 = NodeGroupLogger()
log2 = NodeGroupLogger()
loggers = [:(y == 1.0) << log1, :(x == 0.5) << log2]


setloggers!(model, loggers)


bcs = [
    :(x == 0) << NodeBC(ut = 10.0),
    :(x == 1) << NodeBC(ut = 70.0),
    :(y == 0) << NodeBC(ut = 20.0),
    :(y == 2) << NodeBC(ut = 50.0),

]

tm_solve!(model,end_time = 500000.0, tol = 0.1, nincs = 40)
# Output
save(model, "dom10.vtk")
save(log1, "book3y.dat")
save(log2, "book3x.dat")
