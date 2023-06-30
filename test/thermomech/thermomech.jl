using Amaru
#using Test

# Mesh generation
blocks = [Block([0 0; 1 2], nx = 2, ny = 4, tag = "solids")]

msh = Mesh(blocks)

# Finite element analysis

# Analysis data
k     = 0.0502 # thermal conductivity kW/m/K
rho   = 7.8    # material specific weight Ton/m3
cv    = 486.0  # specific heat (capacity) kJ/Ton/K
E     = 200e6  # kPa
nu    = 0.3
alpha = 1.2e-5 # thermal expansion coefficient  1/K or 1/°C

materials = ["solids" => ElasticSolidThermo(
    E     = E,
    nu    = nu,
    k     = k,
    rho   = rho,
    cv    = cv,
    alpha = alpha,
) ]

model = FEModel(msh, materials)

loggers = [:(y == 1) => NodeGroupLogger("book3.dat")]
setloggers!(model, loggers)

bcs = [
    :(x == 0) => NodeBC(ut = 10.0),
    :(y == 2) => NodeBC(ut = 20.0),
    :(y == 0) => NodeBC(ux = 0, uy = 0),
    :(x == 1) => NodeBC(fx = 100.0),
]
addstage!(model, bcs, tspan=3000000, nincs=10)

tm_solve!(model, tol=0.1)


struct Properties
    rho

end